Build dependencies: 2016-06-14

---

Adios
- adios github version.

OpenFOAM
- develop branch with extras:
  * feature-iotweaks2


For OpenFOAM: fetch copy from develop.openfoam.com (remote: origin)

    git checkout -b myAdiosMix origin/develop

    git merge origin/feature-iotweaks2

    rebuild OpenFOAM with these changes

As a temporary measure (transition), the adios-writer is built
without any lagrangian support.

---
