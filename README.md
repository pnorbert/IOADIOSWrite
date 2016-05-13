## Preliminary layout for OpenFOAM data within ADIOS

* This is a work-in-progress (proof-of-concept) for reading/writing
  OpenFOAM data with ADIOS.
* For development purposes, a function object is used for the implementation framework.
* Restart is with a 'fast-forward' principle: OpenFOAM is used for the Time=0
  construction of meshes and fields, fields are overwritten with their
  values read from the adios file.
* Supports multiple regions, multiple clouds per region.
* Restarting with mesh changing is not yet implemented.
* Restarting with clouds needs more work.

---

2016-05-13

---

* Storage directory:  "adiosData/"

* Each "bp" (binary packed) adios file is restricted to a ***single***
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
All entries are considered mandatory.


| type    | name                        | example
|---------|-----------------------------|--------
| _any_   | /constant/...               | _reserved_
| _any_   | /system/...                 | _reserved_
| string  | /openfoam/version           | "plus-e40d8870f95a"
| string  | /openfoam/platform          | "linux64Gcc"
| int     | /openfoam/label-size        | 32
| string  | /openfoam/precision         | "double"
| int     | /openfoam/nProcs            | 4
| int     | /openfoam/nRegions          | 2
| string[]| /openfoam/regions           | {"region0", "solid"}
| string[]| /openfoam/moving            | {"solid"}
| string[]| /openfoam/topo-change       | {"region0"}


The number of regions can either be obtained directly from the
`/openfoam/nRegions` attribute, or aternatively from the list length
of the `/openfoam/regions` attribute. No particular sort order is
specified for the region names.

### General Variables

Coordinated data for all mesh, fields and cloud information.

| type    | name                        | comment
|---------|-----------------------------|-------------
| _any_   | /constant/...               | _reserved_
| _any_   | /system/...                 | _reserved_


### Time Attributes

Since each file is restricted to a ***single*** time-step/iteration,
which inherently identical across all processes, the time management
values are tracked as attributes rather than as variables. This
additionally simplifies later post-processing, since the attribute
values are available directly from the meta-data.


| type    | name            | comment
|---------|-----------------|-------------
| int     | /time/index     | iteration value
| double  | /time/value     | time-value
| double  | /time/deltaT    | current time-step value
| double  | /time/deltaT0   | previous time-step value

The attribute `/time/value` can be considered to be a global time value
for all fields and clouds contained within the file.


### Region Attributes

Basic information that applies to mesh and field information.

| type    | name                        | comment
|---------|-----------------------------|--------
| int     | \<regionName\>/nPatches     | number of patches
| string[]| \<regionName\>/patch-names  | patch names
| string[]| \<regionName\>/patch-types  | patch types

#### Example

| type    | name                        | example
|---------|-----------------------------|--------
| int     | region0/nPatches            | 3
| string[]| region0/patch-names         | {"movingWall", "fixedWalls", "frontAndBack"}
| string[]| region0/patch-types         | {"wall", "wall", "empty" }
|         | solid/...                   |


### Meshes

These variables are only available when the *updateMesh* attribute is
also true.

#### Global Array Variables


| type    | name                                | comment
|---------|-------------------------------------|-------------
| int     | \<regionName\>/polyMesh/nPoints     | per processor
| int     | \<regionName\>/polyMesh/nCells      | per processor
| int     | \<regionName\>/polyMesh/nFaces      | per processor
| int     | \<regionName\>/polyMesh/nInternalFaces  | per processor


#### Localized Variables

| type    | name                                  | comment
|---------|---------------------------------------|-------------
| double  | \<regionName\>/polyMesh/points        | {nPoints x 3}
| int     | \<regionName\>/polyMesh/faces/indices | indices for compact faceList
| int     | \<regionName\>/polyMesh/faces/content | content for compact faceList
| int     | \<regionName\>/polyMesh/owner         | face owners
| int     | \<regionName\>/polyMesh/neighbour     | face neighbours
| byte    | \<regionName\>/polyMesh/boundary      | boundary mesh information


### Fields


#### Localized Variables

| type    | name                                | comment
|---------|-------------------------------------|-------------
| byte    | \<regionName\>/field/p              |
| byte    | \<regionName\>/field/U              |
|         | ...


#### Attributes

| type    | name                                | example
|---------|-------------------------------------|-------------
| string  | \<regionName\>/field/p/class        | "volScalarField"
| string[]| \<regionName\>/field/p/patch-types  | { "zeroGradient", "zeroGradient", "empty" }
| string  | \<regionName\>/field/U/class        | "volVectorField"
| string[]| \<regionName\>/field/U/patch-types  | {"fixedValue", "fixedValue", "empty"}
|         | ...


### Clouds

These entries are only available when the region also has any active
cloud information. A missing value of `nClouds` is equivalent to `nClouds=0`.

For efficiency, the parcel contents are stored directly as a binary
*blob*. The meaning of the binary content can be reconstructed
from its list of names and types.

#### Region Cloud Attributes

| type    | name                        | comment
|---------|-----------------------------|-------------
| int     | \<regionName\>/nClouds      | number of associated clouds in region
| string[]| \<regionName\>/clouds       | {"name0", "name1", ...}


#### Cloud Variables (Globally-addressable)

| type    | name                           | comment
|---------|--------------------------------|-------------
| byte    | \<regionName\>/cloud/\<cloudName\> | {nParticle x Blob-size}
| int     | \<regionName\>/cloud/\<cloudName\>/nParticle | per processor


#### Cloud Attributes

| type    | name                           | example
|---------|--------------------------------|-------------
| int     | \<regionName\>/cloud/\<cloudName\>/nParticle | total count (sum of corresponding variable)
| string  | \<regionName\>/cloud/\<cloudName\>/class  | "Cloud\<basicKinematicCollidingParcel\>"
| int     | \<regionName\>/cloud/\<cloudName\>/size   | 164 (bytes)
| string[]| \<regionName\>/cloud/\<cloudName\>/names  | {"position", "cellI", "faceI", "stepFraction", "tetFaceI", "tetPtI", "origProc", "origId", "active", "typeId", "nParticle", "d", "dTarget", "U", "rho", "age", "tTurb", "UTurb"}
| string[]| \<regionName\>/cloud/\<cloudName\>/types  | {"vector", "label", "label", "scalar", "label", "label", "label", "label", "label", "label", "scalar", "scalar", "scalar", "vector", "scalar", "scalar", "scalar", "vector"}
| int[]   | \<regionName\>/cloud/\<cloudName\>/offset | {1, 25, 29, 33, 41, 45, 49, 53, 59, 63, 67, 75, 83, 91, 115, 123, 131, 139}
| int[]   | \<regionName\>/cloud/\<cloudName\>/byte-size | {24, 4, 4, 8, 4, 4, 4, 4, 4, 4, 8, 8, 8, 24, 8, 8, 8, 24}

---
