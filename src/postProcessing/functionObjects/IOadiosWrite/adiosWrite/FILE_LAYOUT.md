## Preliminary layout for OpenFOAM data within ADIOS

2016-05-06

---

* Storage directory:  "adiosData/"

* Each "bp" (binary packed) adios file is restricted to a single
  time-step/iteration. This makes for simple and efficient handling.

* In rare cases are values stored in global arrays with offsets.
  - Single values per-processor (eg, time index).
  - Cloud parcel data (currently).

  In all other cases, data are stored with their local dimensions only.

* Data stored as "unsigned byte" generally represent binary content
  that has been serialized via the OpenFOAM `Ostream` and are targeted
  for use by an `Istream` when reading. This binary content should be
  largely identical to a normal OpenFOAM field-file, but without the
  file-header.

---

### General Attributes

These attributes provide assistance when reading the data files.


| type    | name                        | example
|---------|-----------------------------|--------
| _any_   | /constant/...               | _reserved_
| _any_   | /system/...                 | _reserved_
| string  | /openfoam/version           | "plus-e40d8870f95a"
| string  | /openfoam/platform          | "linux64Gcc"
| int     | /openfoam/label-size        | 32
| string  | /openfoam/precision         | "double"
| int     | /openfoam/updateMesh        | 1 (bool, but may change in the future to include indexing etc)
| int     | /openfoam/nProcs            | 4
| int     | /openfoam/nRegions          | 2
| string[]| /openfoam/regions           | {"region0", "solid"}


### General Variables

Coordinated data for all mesh, fields and cloud information.

| type    | name                        | comment
|---------|-----------------------------|-------------
| _any_   | /constant/...               | _reserved_
| _any_   | /system/...                 | _reserved_
| int     | /time/index                 | iteration
| double  | /time/value                 | time-value
| double  | /time/deltaT                |
| double  | /time/deltaT0               |



### Region Attributes

Basic information that applies to mesh and field information.

| type      | name                          | example
|-----------|-------------------------------|--------
| string    | region0/name                  | "region0"
| integer   | region0/npatches              | 3
| string    | region0/patch0/name           | "movingWall"
| string    | region0/patch0/type           | "wall"
| string    | region0/patch1/name           | "fixedWalls"
| string    | region0/patch1/type           | "wall"
| string    | region0/patch2/name           | "frontAndBack"
| string    | region0/patch2/type           | "empty"
|           | region1/...                   |


### Meshes

These variables are only available when the *updateMesh* attribute is
also true.

#### Global Array Variables


| type      | name                          | comment
|-----------|-------------------------------|-------------
| integer   | region0/polyMesh/nPoints      | per processor
| integer   | region0/polyMesh/nCells       | per processor
| integer   | region0/polyMesh/nFaces       | per processor
| integer   | region0/polyMesh/nInternalFaces | per processor
|           | region1/...                   |


#### Localized Variables

| type      | name                           | comment
|-----------|--------------------------------|-------------
| double    | region0/polyMesh/points        | {nPoints x 3}
| integer   | region0/polyMesh/faces/indices | indices for compact faceList
| integer   | region0/polyMesh/faces/content | content for compact faceList
| integer   | region0/polyMesh/owner         | face owners
| integer   | region0/polyMesh/neighbour     | face neighbours
| byte      | region0/polyMesh/boundary      | boundary mesh information
|           | region1/...                    |


### Fields


#### Localized Variables

| type      | name                           | comment
|-----------|--------------------------------|-------------
| byte      | region0/field/p                |
| byte      | region0/field/U                |
|           | ...


#### Attributes

| type      | name                           | example
|-----------|--------------------------------|-------------
| string    | region0/field/p/class          | "volScalarField"
| string    | region0/field/p/patch0/type    | "zeroGradient"
| string    | region0/field/p/patch1/type    | "zeroGradient"
| string    | region0/field/p/patch2/type    | "empty"
| string    | region0/field/U/class          | "volVectorField"
| string    | region0/field/U/patch0/type    | "fixedValue"
| string    | region0/field/U/patch1/type    | "fixedValue"
| string    | region0/field/U/patch2/type    | "empty"
|           | ...


### Clouds

For efficiency, the parcel contents are stored directly as a binary
*blob*. The meaning of the binary content can be reconstructed
from its list of names and types.

#### Global Array Variables

| type      | name                           | comment
|-----------|--------------------------------|-------------
| integer   | region0/cloud0/nParticle       | per processor
| byte      | region0/cloud0/\_\_blob\_\_    | {nParticle x Blob-size}

#### Cloud Attributes

| type      | name                           | example
|-----------|--------------------------------|-------------
| integer   | region0/nClouds                | 1
| string    | region0/cloud0/class           | "kinematicCloud"
| string    | region0/cloud0/name            | "kinematicCloud"
| integer   | region0/cloud0/nParticle       | 3000 (total count)
| char*[]   | region0/cloud0/\_\_blob\_\_/names  | {"position", "cellI", "faceI", "stepFraction", "tetFaceI", "tetPtI", "origProc", "origId", "active", "typeId", "nParticle", "d", "dTarget", "U", "rho", "age", "tTurb", "UTurb"}
| char*[]   | region0/cloud0/\_\_blob\_\_/types  | {"vector", "label", "label", "scalar", "label", "label", "label", "label", "label", "label", "scalar", "scalar", "scalar", "vector", "scalar", "scalar", "scalar", "vector"}
| integer   | region0/cloud0/\_\_blob\_\_/offset | {0, 24, 28, 32, 40, 44, 48, 52, 56, 60, 64, 72, 80, 88, 112, 120, 128, 136}
| integer   | region0/cloud0/\_\_blob\_\_/byte-size | {24, 4, 4, 8, 4, 4, 4, 4, 4, 4, 8, 8, 8, 24, 8, 8, 8, 24}


---
