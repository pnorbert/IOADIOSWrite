/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |               2015 Norbert Podhorszki
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "adiosWrite.H"

#include "basicKinematicCloud.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::adiosWrite::cloudDefine(label regionID)
{
    regionInfo& r = regions_[regionID];
    Info<< "  adiosWrite::cloudDefine: region " << regionID << " " << r.name_ << ": " << endl;

    const fvMesh& mesh = time_.lookupObject<fvMesh>(regions_[regionID].name_);

    forAll(r.cloudNames_, cloudI)
    {
        Info<< "    cloud: " << r.cloudNames_[cloudI] << endl;

        const kinematicCloud& cloud =
            mesh.lookupObject<kinematicCloud>(r.cloudNames_[cloudI]);
        const basicKinematicCloud& q =
            static_cast<const basicKinematicCloud&>(cloud);

        // Number of particles on this process
        label myParticles = q.size();

        // Find the number of particles on each process
        r.nParticles_[Pstream::myProcNo()] = myParticles;
        Pstream::gatherList(r.nParticles_);
        Pstream::scatterList(r.nParticles_);

        // Sum total number of particles on all processes
        r.nTotalParticles_ = sum(r.nParticles_);

        // If the cloud contains no particles, jump to the next cloud
        if (r.nTotalParticles_ == 0)
        {
            Info<< "    " << r.cloudNames_[cloudI] <<": No particles in cloud. "
                << "Skipping definition." <<endl;
            continue;
        }

        // Calculate offset values
        List<label> offsets(Pstream::nProcs());
        offsets[0] = 0;
        for (label proc=1; proc < offsets.size(); proc++)
        {
            offsets[proc] = offsets[proc-1] + r.nParticles_[proc-1];
        }

        //List<word> intsets(4) = {"origProc", "origID", "cell", "currProc"};

        // Define a variable for dataset name
        char varPath[256];
        char gdimstr[16]; // 1D global array of particles from all processes
        char ldimstr[16]; // this process' size
        char offsstr[16]; // this process' offset in global array

        sprintf (gdimstr, "%d", r.nTotalParticles_);
        sprintf (ldimstr, "%d", myParticles);
        sprintf (offsstr, "%d", offsets[Pstream::myProcNo()]);

        // Define all possible output variables, not necessary to write all of them later
        sprintf (varPath, "region%d/clouds/%s", regionID, r.cloudNames_[cloudI].c_str());

        adios_define_var (groupID_, "nParticlesPerProc", varPath, adios_integer, NULL, NULL, NULL);
        adios_define_var (groupID_, "nTotalParticles", varPath, adios_integer, NULL, NULL, NULL);
        outputSize_ += 2*sizeof(int);

        // integer info datasets
        adios_define_var (groupID_, "origProc", varPath, adios_integer, ldimstr, gdimstr, offsstr);
        adios_define_var (groupID_, "origId",   varPath, adios_integer, ldimstr, gdimstr, offsstr);
        adios_define_var (groupID_, "cell",     varPath, adios_integer, ldimstr, gdimstr, offsstr);
        adios_define_var (groupID_, "currProc", varPath, adios_integer, ldimstr, gdimstr, offsstr);
        outputSize_ += 4 * myParticles * sizeof(int);

        // scalar datasets
        adios_define_var (groupID_, "rho", varPath, ADIOS_SCALAR, ldimstr, gdimstr, offsstr);
        adios_define_var (groupID_, "d",   varPath, ADIOS_SCALAR, ldimstr, gdimstr, offsstr);
        adios_define_var (groupID_, "age", varPath, ADIOS_SCALAR, ldimstr, gdimstr, offsstr);
        outputSize_ += 3 * myParticles * sizeof(ioScalar);

        // vector datasets
        sprintf (gdimstr, "%d,3", r.nTotalParticles_);
        sprintf (ldimstr, "%d,3", myParticles);
        sprintf (offsstr, "%d,0", offsets[Pstream::myProcNo()]);
        adios_define_var (groupID_, "position", varPath, ADIOS_SCALAR, ldimstr, gdimstr, offsstr);
        adios_define_var (groupID_, "U",        varPath, ADIOS_SCALAR, ldimstr, gdimstr, offsstr);
        adios_define_var (groupID_, "Us",       varPath, ADIOS_SCALAR, ldimstr, gdimstr, offsstr);
        outputSize_ += 3 * myParticles * sizeof(ioScalar) * 3;
    }
}


void Foam::adiosWrite::cloudWrite(label regionID)
{
    const regionInfo& r = regions_[regionID];
    const fvMesh& mesh = time_.lookupObject<fvMesh>(r.name_);
    Info<< "  adiosWrite::cloudWrite: region " << regionID << " " << r.name_ << ": " << endl;

    forAll(r.cloudNames_, cloudI)
    {
        Info<< "    cloud: " << r.cloudNames_[cloudI] << endl;

        // If the cloud contains no particles, jump to the next cloud
        if (r.nTotalParticles_ == 0)
        {
            continue;
        }

        const kinematicCloud& cloud =
            mesh.lookupObject<kinematicCloud>(r.cloudNames_[cloudI]);
        const basicKinematicCloud& q =
            static_cast<const basicKinematicCloud&>(cloud);

        label myParticles = r.nParticles_[Pstream::myProcNo()];

        char varPath[256], datasetName[256];
        sprintf (varPath, "region%d/clouds/%s", regionID, r.cloudNames_[cloudI].c_str());

        // Allocate memory for 1-comp. dataset of type 'integer' for adios_integer writes
        int particleLabel[myParticles];

        sprintf (datasetName, "%s/nParticlesPerProc", varPath);
        int nparts = myParticles;
        int ntotalparts = r.nTotalParticles_;
        adios_write (fileID_, datasetName, &nparts);
        sprintf (datasetName, "%s/nTotalParticles", varPath);
        adios_write (fileID_, datasetName, &ntotalparts);

        // Write original processor ID
        if (findStrings(r.cloudAttribs_, "origProc"))
        {
            Info<< "      dataset origProc " << endl;
            label i = 0;
            forAllConstIter(basicKinematicCloud, q, pIter)
            {
                particleLabel[i++] = pIter().origProc();
            }
            sprintf (datasetName, "%s/origProc", varPath);
            adios_write (fileID_, datasetName, particleLabel);
        }

        // Write original ID
        if (findStrings(r.cloudAttribs_, "origId"))
        {
            Info<< "      dataset origId " << endl;
            label i = 0;
            forAllConstIter(basicKinematicCloud, q, pIter)
            {
                particleLabel[i++] = pIter().origId();
            }
            sprintf (datasetName, "%s/origId", varPath);
            adios_write (fileID_, datasetName, particleLabel);
        }

        // Write cell number
        if (findStrings(r.cloudAttribs_, "cell"))
        {
            Info<< "      dataset cell " << endl;
            label i = 0;
            forAllConstIter(basicKinematicCloud, q, pIter)
            {
                particleLabel[i++] = pIter().cell();
            }
            sprintf (datasetName, "%s/cell", varPath);
            adios_write (fileID_, datasetName, particleLabel);
        }

        // Write current process ID
        if (findStrings(r.cloudAttribs_, "currProc"))
        {
            Info<< "      dataset currProc " << endl;
            label i = 0;
            forAllConstIter(basicKinematicCloud, q, pIter)
            {
                particleLabel[i++] = Pstream::myProcNo();
            }
            sprintf (datasetName, "%s/currProc", varPath);
            adios_write (fileID_, datasetName, particleLabel);
        }

        // Allocate memory for 1-comp. dataset of type 'ioScalar'
        ioScalar* particleScalar1;
        particleScalar1 = new ioScalar[myParticles];

        // Write density rho
        if (findStrings(r.cloudAttribs_, "rho"))
        {
            Info<< "      dataset rho " << endl;
            label i = 0;
            forAllConstIter(basicKinematicCloud, q, pIter)
            {
                particleScalar1[i++] = pIter().rho();
            }
            sprintf (datasetName, "%s/rho", varPath);
            adios_write (fileID_, datasetName, particleScalar1);
        }

        // Write diameter d
        if (findStrings(r.cloudAttribs_, "d"))
        {
            Info<< "      dataset d " << endl;
            label i = 0;
            forAllConstIter(basicKinematicCloud, q, pIter)
            {
                particleScalar1[i++] = pIter().d();
            }
            sprintf (datasetName, "%s/d", varPath);
            adios_write (fileID_, datasetName, particleScalar1);
        }

        // Write age
        if (findStrings(r.cloudAttribs_, "age"))
        {
            label i = 0;
            forAllConstIter(basicKinematicCloud, q, pIter)
            {
                particleScalar1[i++] = pIter().age();
            }
            sprintf (datasetName, "%s/age", varPath);
            adios_write (fileID_, datasetName, particleScalar1);
        }

        // Free memory for 1-comp. dataset of type 'ioScalar'
        delete [] particleScalar1;


        // Allocate memory for 3-comp. dataset of type 'ioScalar'
        ioScalar* particleScalar3;
        particleScalar3 = new ioScalar[myParticles*3];

        // Write position
        if (findStrings(r.cloudAttribs_, "position"))
        {
            label i = 0;
            forAllConstIter(basicKinematicCloud, q, pIter)
            {
                particleScalar3[3*i+0] = pIter().position().x();
                particleScalar3[3*i+1] = pIter().position().y();
                particleScalar3[3*i+2] = pIter().position().z();
                i++;
            }
            sprintf (datasetName, "%s/position", varPath);
            adios_write (fileID_, datasetName, particleScalar3);
        }

        // Write velocity U
        if (findStrings(r.cloudAttribs_, "U"))
        {
            label i = 0;
            forAllConstIter(basicKinematicCloud, q, pIter)
            {
                particleScalar3[3*i+0] = pIter().U().x();
                particleScalar3[3*i+1] = pIter().U().y();
                particleScalar3[3*i+2] = pIter().U().z();
                i++;
            }
            sprintf (datasetName, "%s/U", varPath);
            adios_write (fileID_, datasetName, particleScalar3);
        }

        // Write slip velocity Us = U - Uc
        if (findStrings(r.cloudAttribs_, "Us"))
        {
            label i = 0;
            forAllConstIter(basicKinematicCloud, q, pIter)
            {
                particleScalar3[3*i+0] = pIter().U().x() - pIter().Uc().x();
                particleScalar3[3*i+1] = pIter().U().y() - pIter().Uc().y();
                particleScalar3[3*i+2] = pIter().U().z() - pIter().Uc().z();
                i++;
            }
            sprintf (datasetName, "%s/Us", varPath);
            adios_write (fileID_, datasetName, particleScalar3);
        }

        // Free memory for 3-comp. dataset of type 'ioScalar'
        delete [] particleScalar3;
    }
}


// ************************************************************************* //
