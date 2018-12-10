# General Conjugate Gradient Eigensolver (GCGE)

## OVERVIEW

| Document Name | Description |
| ------------- | ----------- |
| src           |             |
| app           |             |
| test          |             |
| blaslapack    |             |


## INSTALLATION DETAILS

We uses a simplified `Makefile` to compile and install the package. 

User should set up compiler names (FORTRAN and C) and compilation options 

in the file `config/make.inc`. Type `make help` for more details.

`$ vi config/make.inc`

provide correct names and options

`$ make libs`

`$ make install`

## RUNNING DEMOs

## STRUCTURE

> src/
>
> > Makefile
> > gcge.h
> > gcge_cg.c
> > gcge_cg.h
> > gcge_eigsol.c
> > gcge_eigsol.h
> > gcge_matvec.c
> > gcge_matvec.h
> > gcge_ops.c
> > gcge_ops.h
> > gcge_orthogonal.c
> > gcge_orthogonal.h
> > gcge_para.c
> > gcge_para.h
> > gcge_rayleighritz.c
> > gcge_rayleighritz.h
> > gcge_solver.c
> > gcge_solver.h
> > gcge_statistics.c
> > gcge_statistics.h
> > gcge_type.h
> > gcge_utilities.c
> > gcge_utilities.h
> > gcge_workspace.c
> > gcge_workspace.h
> > gcge_xpw.c
> > gcge_xpw.h

> app/
>
> > csr
> > > Makefile
> > > gcge_app_csr.c
> > > gcge_app_csr.h
>
> > hypre
> > > Makefile
> > > gcge_app_hypre.c
> > > gcge_app_hypre.h
>
> > petsc
> > > Makefile
> > > gcge_app_petsc.c
> > > gcge_app_petsc.h
>
> > phg
> > > Makefile
> > > gcge_app_phg.c
> > > gcge_app_phg.h
>
> > slepc
> > > Makefile
> > > gcge_app_slepc.c
> > > gcge_app_slepc.h
