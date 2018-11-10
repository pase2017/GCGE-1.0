/*
 * =====================================================================================
 *
 *       Filename:  test_solver.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2018年09月24日 09时57分13秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

//#include "memwatch.h"

#include "gcge.h"
#include "gcge_app_phg.h"

#include "phg.h"

#if (PHG_VERSION_MAJOR <= 0 && PHG_VERSION_MINOR < 9)
# undef ELEMENT
typedef SIMPLEX ELEMENT;
#endif

/* analytic solution of the first eigen vector in the case of the unit cube */
static int nev_plus, (*modes)[3], *mode;

static void
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = Sin(mode[0] * M_PI * x) *
	     Sin(mode[1] * M_PI * y) *
	     Sin(mode[2] * M_PI * z);
} 

FLOAT
compute_error(DOF *u_h, int which)
{
    /* current range is [start, start + n) */
    static int start = -1, n = 0, *pvt = NULL;
    static FLOAT *A = NULL, *b;
    static DOF **v = NULL;
    DOF *u;
    INT i, k;

    if (u_h == NULL) {
	/* free buffers */
	if (v != NULL) {
	    for (i = 0; i < n; i++)
		phgDofFree(v + i);
	    FreeAtExit(A);
	    FreeAtExit(pvt);
	    FreeAtExit(v);
	    A = NULL;
	    pvt = NULL;
	    v = NULL;
	}
	start = -1;
	n = 0;
	return 0.;
    }

    mode = (int *)(modes + which);
    i = mode[0] * mode[0] + mode[1] * mode[1] + mode[2] * mode[2];
    for (k = which; k >= 0; k--) {
	mode = (int *)(modes + k);
	if (i != mode[0] * mode[0] + mode[1] * mode[1] + mode[2] * mode[2])
	    break;
    }
    if ((++k) != start) {
	if (start >= 0 && n > 0) {
	    assert(v != NULL);
	    for (start = 0; start < n; start++)
		phgDofFree(v + start);
	    phgFree(A);
	    phgFree(pvt);
	    phgFree(v);
	}
	start = k;
	for (n = which + 1; n < nev_plus; n++) {
	    mode = (int *)(modes + n);
	    if (i != mode[0] * mode[0] + mode[1] * mode[1] + mode[2] * mode[2])
		break;
	}
	n -= start;
	A = phgAlloc((n + 1) * n * sizeof(*A));
	b = A + n * n;
	pvt = phgAlloc(n * sizeof(*pvt));
	v = phgAlloc(n * sizeof(*v));
	for (i = 0; i < n; i++) {
	    mode = (int *)(modes + start + i);
	    v[i] = phgDofNew(u_h->g, u_h->type, 1, "v", func_u);
	}
	for (i = 0; i < n; i++) {
	    A[i * n + i] = phgDofDotL2Vec(v[i], v[i]);
	    for (k = i + 1; k < n; k++)
		A[i * n + k] = A[k * n + i] = phgDofDotL2Vec(v[i], v[k]);
	}
	phgSolverDenseLU(n, A, pvt);
    }

    for (i = 0; i < n; i++)
	b[i] = phgDofDotL2Vec(v[i], u_h);
    phgSolverDenseSV(n, A, pvt, 1, b);
    u = NULL;
    for (i = 0; i < n; i++)
	phgDofAXPY(b[i], v[i], &u);
    phgDofAXPY(-1.0, u_h, &u);
    b[0] = phgDofNormL2Vec(u) / phgDofNormL2Vec(u_h);
    phgDofFree(&u);

    return b[0];
}

static void
build_matrices(MAT *matA, MAT *matB, DOF *u_h)
{
    int N = u_h->type->nbas;	/* number of basis functions in an element */
    int i, j;
    GRID *g = u_h->g;
    ELEMENT *e;
    FLOAT A[N][N], B[N][N], vol;
    static FLOAT *B0 = NULL;
    INT I[N];

    if (B0 == NULL && g->nroot > 0) {
	/* (\int \phi_j\cdot\phi_i)/vol is independent of element */
	FreeAtExit(B0);
	B0 = phgAlloc(N * N * sizeof(*B0));
	e = g->roots;
	vol = phgGeomGetVolume(g, e);
	for (i = 0; i < N; i++)
	    for (j = 0; j <= i; j++)
		B0[i * N + j] = B0[i + j * N] =
		    phgQuadBasDotBas(e, u_h, j, u_h, i, QUAD_DEFAULT) / vol;
    }

    assert(u_h->dim == 1);
    ForAllElements(g, e) {
	vol = phgGeomGetVolume(g, e);
	for (i = 0; i < N; i++) {
	    I[i] = phgMapE2L(matA->cmap, 0, e, i);
	    for (j = 0; j <= i; j++) {
		/* \int \grad\phi_j\cdot\grad\phi_i */
		A[j][i] = A[i][j] =
		    phgQuadGradBasDotGradBas(e, u_h, j, u_h, i, QUAD_DEFAULT);
		/* \int \phi_j\cdot\phi_i */
		B[j][i] = B[i][j] = B0[i * N + j] * vol;
	    }
	}

	/* loop on basis functions */
	for (i = 0; i < N; i++) {
	    if (phgDofDirichletBC(u_h, e, i, NULL, NULL, NULL, DOF_PROJ_NONE))
		continue;
	    phgMatAddEntries(matA, 1, I + i, N, I, A[i]); 
	    phgMatAddEntries(matB, 1, I + i, N, I, B[i]); 
	}
    }
}

static FLOAT
estimate_error(FLOAT lambda, DOF *u_h, DOF *error)
/* compute H1 error indicator */
{
    GRID *g = u_h->g;
    ELEMENT *e;
    DOF *grad_u, *jump, *residual, *tmp;
    FLOAT PEoo = 0.0;

    grad_u = phgDofGradient(u_h, NULL, NULL, NULL);
    tmp = phgDofDivergence(grad_u, NULL, NULL, NULL);
    residual = phgDofGetSameOrderDG(u_h, -1, NULL);
    phgDofCopy(u_h, &residual, NULL, NULL);
    phgDofAXPBY(1.0, tmp, lambda, &residual);
    phgDofFree(&tmp);
    jump = phgQuadFaceJump(grad_u, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
    ForAllElements(g, e) {
	int i;
	FLOAT eta, h;
	FLOAT diam = phgGeomGetDiameter(g, e);
	eta = 0.0;
	/* for each face F compute [grad_u \cdot n] */
	for (i = 0; i < NFace; i++) {
	    if (e->bound_type[i] & (DIRICHLET | NEUMANN))
		continue;	/* boundary face */
	    h = phgGeomGetFaceDiameter(g, e, i);
	    eta += *DofFaceData(jump, e->faces[i]) * h;
	}
	eta = eta * .5 + diam*diam*phgQuadDofDotDof(e, residual, residual, -1);
	if (*DofElementData(error, e->index) < eta)
	    *DofElementData(error, e->index) = eta;
	if (PEoo < eta)
	    PEoo = eta;
    }
    phgDofFree(&jump);
    phgDofFree(&residual);
    phgDofFree(&grad_u);

    return Sqrt(PEoo);
}

static int
bc_map(int bctype)
{
    return DIRICHLET;	/* set Dirichlet BC on all boundaries */
}

/* Note: make USER_CFLAGS="-DMATRIX_FREE_TEST=1" eigen */
#ifndef MATRIX_FREE_TEST
# define MATRIX_FREE_TEST 0
#endif	/* !defined(MATRIX_FREE_TEST) */

#if MATRIX_FREE_TEST
/* callback function for matrix A */
static int
mat_vec(MAT_OP op, MAT *mat, VEC *x, VEC *y)
/* computes y = op(A) * x */
{
    phgMatVec(op, 1.0, ((MAT **)mat->mv_data)[0], x, 0.0, &y);
    return 0;
}
#endif	/* MATRIX_FREE_TEST */

int
main(int argc, char *argv[])
{
    static char *fn = "../data/cube4.dat";
    /*static char *fn = "../test/model.mesh";*/
    static INT mem_max = 400;
    size_t mem, mem_peak;
    int i, j, k, n, nit;
    INT nev = 10;
    INT pre_refines = 2;
    GRID *g;
    DOF **u_h, *error;
    MAP *map;
    MAT *A, *B;
    FLOAT tol = 1e-3, tau = 0.0, PEoo, thres;
    FLOAT *evals;
    double wtime;
#if MATRIX_FREE_TEST
    MAT *C;

    phgOptionsPreset("-arpack_solver pcg");
#endif	/* MATRIX_FREE_TEST */

    phgOptionsPreset("-dof_type P4");

    phgOptionsRegisterFilename("-mesh_file", "Mesh filename", &fn);
    phgOptionsRegisterInt("-pre_refines", "Pre-refines", &pre_refines);
    phgOptionsRegisterInt("-nev", "Number of eigenvalues", &nev);
    phgOptionsRegisterFloat("-tol", "Convergence tolerance", &tol);
    phgOptionsRegisterFloat("-tau", "The shift", &tau);
    phgOptionsRegisterInt("-mem_max", "Maximum memory (MB)", &mem_max);

    phgInit(&argc, &argv);

    if (DOF_DEFAULT->mass_lumping == NULL)
	phgPrintf("Order of FE bases: %d\n", DOF_DEFAULT->order);
    else
	phgPrintf("Order of FE bases: %d\n", DOF_DEFAULT->mass_lumping->order0);

    g = phgNewGrid(-1);
    phgImportSetBdryMapFunc(bc_map);
    /*phgSetPeriodicity(g, X_MASK | Y_MASK | Z_MASK);*/
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);
    /* pre-refinement */
    phgRefineAllElements(g, pre_refines);

    u_h = phgAlloc(nev * sizeof(*u_h));
    evals = phgCalloc(nev, sizeof(*evals));
    for (i = 0; i < nev; i++) {
	u_h[i] = phgDofNew(g, DOF_DEFAULT, 1, "u_h", DofInterpolation);
#if 0
	/* All-Neumann BC */
	phgDofSetDirichletBoundaryMask(u_h[i], 0);
#else
	/* All-Dirichlet BC */
	phgDofSetDirichletBoundaryMask(u_h[i], BDRY_MASK);
#endif
	/* set random initial values for the eigenvectors */
	phgDofRandomize(u_h[i], i == 0 ? 123 : 0);
    }
    error = phgDofNew(g, DOF_P0, 1, "error indicator", DofNoAction);

    nev_plus = nev + 1000;
    modes = phgAlloc(nev_plus * sizeof(modes[0]));
    for (n = 0, nit = 3; n < nev_plus; nit++) {
	/* find i, j, k such that i^2 + j^2 + k^2 == nit */
	for (i = 1; n < nev_plus && i * i <= nit; i++) {
	    for (j = 1; n < nev_plus && j * j <= nit - i * i; j++) {
		for (k = 1; k * k < nit - i * i - j * j; k++);
		if (i * i + j * j + k * k != nit)
		    continue;
		modes[n][0] = i;
		modes[n][1] = j;
		modes[n][2] = k;
		n++;
	    }
	}
    }

#if 0 * USE_BLOPEX
const char *pc1 = phgOptionsGetKeyword("-blopex_pc_solver1");
/* disable local PC for the first loop */
phgOptionsSetKeyword("-blopex_pc_solver1", "none");
#else
const char *pc1 = NULL;
#endif
    while (TRUE) {
	phgPrintf("\n");
	if (phgBalanceGrid(g, 1.2, -1, NULL, 0.))
	    phgPrintf("Repartition mesh\n");
	phgPrintf("%"dFMT" DOF, %"dFMT
		  " elements, %d submesh%s, load imbalance: %lg\n",
			DofGetDataCountGlobal(u_h[0]), g->nleaf_global,
			g->nprocs, g->nprocs > 1 ? "es" : "", (double)g->lif);
	wtime = phgGetTime(NULL);
	map = phgMapCreate(u_h[0], NULL);
#if 0
	A = phgMapCreateMat(map, map);
	B = phgMapCreateMat(map, map);
	build_matrices(A, B, u_h[0]);
	phgMatRemoveBoundaryEntries(A);
	phgMatRemoveBoundaryEntries(B);
#else
	{
	    MAP *m = phgMapRemoveBoundaryEntries(map);
	    A = phgMapCreateMat(m, m);
	    B = phgMapCreateMat(m, m);
	    phgMapDestroy(&m);
	}
	build_matrices(A, B, u_h[0]);
#endif
	phgPrintf("Build matrices, matrix size: %"dFMT", wtime: %0.2lfs\n",
			A->rmap->nglobal, phgGetTime(NULL) - wtime);
	wtime = phgGetTime(NULL);
#if MATRIX_FREE_TEST
	C = phgMapCreateMatrixFreeMat(A->rmap, A->cmap, mat_vec, A, NULL);
	n = phgDofEigenSolve(C, B, nev, EIGEN_CLOSEST, tau, &nit, evals,
			     map, u_h, NULL);
	phgMatDestroy(&C);
#else	/* MATRIX_FREE_TEST */



   //创建一个特征值solver实例
   GCGE_SOLVER *phg_solver = GCGE_PHG_Solver_Init(A, B, nev, argc,  argv);   
   //一些参数的设置
   phg_solver->para->ev_tol = tol*0.01;
   phg_solver->para->dirichlet_boundary = 0;
   //求解特征值问题
   GCGE_SOLVER_Solve(phg_solver);  
   GCGE_SOLVER_Free_All(&phg_solver);


//	n = phgDofEigenSolve(A, B, nev, EIGEN_CLOSEST, tau, &nit, evals,
//			     map, u_h, NULL);





#endif	/* MATRIX_FREE_TEST */
	phgPrintf("%d iterations, converged eigenvalues: %d, wtime: %0.2lfs\n",
			nit, n, phgGetTime(NULL) - wtime);
if (pc1 != NULL && n == nev) phgOptionsSetKeyword("-blopex_pc_solver1", pc1);
	for (i = 0; i < n; i++) {
	    FLOAT err_tau;
	    nit = modes[i][0] * modes[i][0] +
		  modes[i][1] * modes[i][1] +
		  modes[i][2] * modes[i][2];
	    err_tau = Fabs(evals[i] - nit * Pow(M_PI,2)) / evals[i];
	    phgPrintf("  tau[%d]=%0.12e, error(tau)=%0.2e, error(u_h)=%0.2e\n",
			i, (double)evals[i], (double)err_tau,
			(double)compute_error(u_h[i], i));
	}
	compute_error(NULL, 0);		/* reset compute_error() */
	for (; i < nev; i++)
	    phgDofRandomize(u_h[i], 0);
	phgMatDestroy(&B);
	phgMatDestroy(&A);
	phgMapDestroy(&map);

	if (n == 0) {
	    PEoo = 1e10;
	}
	else {
	    PEoo = 0.0;
	    phgDofSetDataByValue(error, 0.0);
	    for (i = 0; i < n; i++) {
		thres = estimate_error(evals[i], u_h[i], error);
		if (PEoo < thres)
		    PEoo = thres;
	    }
	}
	phgPrintf("indicator=%0.3le\n", (double)PEoo);
	mem = phgMemoryUsage(g, &mem_peak);
	phgPrintf("Memory usage: current %0.4lgMB, peak %0.4lgMB\n",
		(double)mem / (1024.0 * 1024.0),
		(double)mem_peak / (1024.0 * 1024.0));
	if (PEoo < tol || mem_peak >= 1024 * (size_t)1024 * mem_max)
	    break;
	if (n == 0) {
	    phgRefineAllElements(g, 1);
	}
	else {
#if 0
	    /*thres = Pow(1e-2,2) / (double)g->nleaf_global;*/
	    thres = Pow(PEoo * 0.5, 2);
	    ForAllElements(g, e)
		if (*DofElementData(error, e->index) > thres)
		    e->mark = 1;
#else
	    phgMarkRefine(MARK_DEFAULT, error, Pow(0.8,2), NULL, 0., 1,
			  Pow(tol, 2) / g->nleaf_global);
#endif
	    phgRefineMarkedElements(g);
	}
    }

#if 1
    phgPrintf("Creating \"%s\".\n", phgExportVTKn(g, "output.vtk", nev, u_h));
#endif

    phgFree(modes);
    phgFree(evals);
    for (i = 0; i < nev; i++) {
	phgDofFree(u_h + i);
    }
    phgDofFree(&error);
    phgFree(u_h);
    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}
