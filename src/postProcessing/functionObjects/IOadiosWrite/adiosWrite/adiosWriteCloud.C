/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |               2015 Norbert Podhorszki
                            |               2016 OpenCFD Ltd.
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

void Foam::adiosWrite::cloudDefine(regionInfo& r)
{
    Info<< "  adiosWrite::cloudDefine: region "
        << r.index_ << " " << r.name_ << ": " << endl;

    const fvMesh& mesh = time_.lookupObject<fvMesh>(r.name_);

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
        char gdimstr[16]; // 1D global array of particles from all processes
        char ldimstr[16]; // this process' size
        char offsstr[16]; // this process' offset in global array

        sprintf(gdimstr, "%d", r.nTotalParticles_);
        sprintf(ldimstr, "%d", myParticles);
        sprintf(offsstr, "%d", offsets[Pstream::myProcNo()]);

        // Define all possible output variables, not necessary to write all of them later
        fileName varPath
        (
            "region" + Foam::name(r.index_)
          / "clouds" / r.cloudNames_[cloudI]
        );

        defineVariable(varPath/"nParticlesPerProc", adios_integer);
        defineVariable(varPath/"nTotalParticles",   adios_integer);

        // integer info datasets
        adios_define_var(groupID_, (varPath/"origProc").c_str(), "", adios_integer, ldimstr, gdimstr, offsstr);
        adios_define_var(groupID_, (varPath/"origId").c_str(),   "", adios_integer, ldimstr, gdimstr, offsstr);
        adios_define_var(groupID_, (varPath/"cell").c_str(),     "", adios_integer, ldimstr, gdimstr, offsstr);
        adios_define_var(groupID_, (varPath/"currProc").c_str(), "", adios_integer, ldimstr, gdimstr, offsstr);
        outputSize_ += 4 * myParticles * sizeof(int);

        // scalar datasets
        adios_define_var(groupID_, (varPath/"rho").c_str(), "", ADIOS_SCALAR, ldimstr, gdimstr, offsstr);
        adios_define_var(groupID_, (varPath/"d").c_str(),   "", ADIOS_SCALAR, ldimstr, gdimstr, offsstr);
        adios_define_var(groupID_, (varPath/"age").c_str(), "", ADIOS_SCALAR, ldimstr, gdimstr, offsstr);
        outputSize_ += 3 * myParticles * sizeof(ioScalar);

        // vector datasets
        sprintf(gdimstr, "%d,3", r.nTotalParticles_);
        sprintf(ldimstr, "%d,3", myParticles);
        sprintf(offsstr, "%d",   offsets[Pstream::myProcNo()]);
        adios_define_var(groupID_, (varPath/"position").c_str(), "", ADIOS_SCALAR, ldimstr, gdimstr, offsstr);
        adios_define_var(groupID_, (varPath/"U").c_str(),        "", ADIOS_SCALAR, ldimstr, gdimstr, offsstr);
        adios_define_var(groupID_, (varPath/"Us").c_str(),       "", ADIOS_SCALAR, ldimstr, gdimstr, offsstr);
        outputSize_ += 3*myParticles*sizeof(ioScalar)*3;
    }
}


void Foam::adiosWrite::cloudWrite(const regionInfo& r)
{
    const fvMesh& mesh = time_.lookupObject<fvMesh>(r.name_);
    Info<< "  adiosWrite::cloudWrite: region "
        << r.index_ << " " << r.name_ << ": " << endl;

    forAll(r.cloudNames_, cloudI)
    {
        Info<< "    cloud: " << r.cloudNames_[cloudI] << endl;
        Info<< "    Properties " <<  basicKinematicCloud::particleType::propertyList() << endl;

        // If the cloud contains no particles, jump to the next cloud
        if (r.nTotalParticles_ == 0)
        {
            continue;
        }

        const kinematicCloud& cloud =
            mesh.lookupObject<kinematicCloud>(r.cloudNames_[cloudI]);
        const basicKinematicCloud& q =
            static_cast<const basicKinematicCloud&>(cloud);

        const label myParticles = r.nParticles_[Pstream::myProcNo()];

        fileName varPath
        (
            "region" + Foam::name(r.index_)
          / "clouds" / r.cloudNames_[cloudI]
        );

        // Allocate storage for 1-comp. dataset of type 'integer'
        int particleLabel[myParticles];
        int nparts = myParticles;
        int ntotalparts = r.nTotalParticles_;

        writeVariable(varPath/"nParticlesPerProc", &nparts);
        writeVariable(varPath/"nTotalParticles",   &ntotalparts);

        // Write original processor ID
        {
            const word what = "origProc";
            if (findStrings(r.cloudAttribs_, what))
            {
                Info<< "      dataset " << what << endl;
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    particleLabel[i++] = pIter().origProc();
                }
                writeVariable(varPath/what, &particleLabel);
            }
        }

        // Write original ID
        {
            const word what = "origId";
            if (findStrings(r.cloudAttribs_, what))
            {
                Info<< "      dataset " << what << endl;
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    particleLabel[i++] = pIter().origId();
                }
                writeVariable(varPath/what, &particleLabel);
            }
        }

        // Write cell number
        {
            const word what = "cell";
            if (findStrings(r.cloudAttribs_, what))
            {
                Info<< "      dataset " << what << endl;
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    particleLabel[i++] = pIter().cell();
                }
                writeVariable(varPath/what, &particleLabel);
            }
        }

        // Write current process ID
        {
            const word what = "currProc";
            if (findStrings(r.cloudAttribs_, what))
            {
                Info<< "      dataset " << what << endl;
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    particleLabel[i++] = Pstream::myProcNo();
                }
                writeVariable(varPath/what, &particleLabel);
            }
        }


        // Allocate storage for 1-comp. dataset of type 'ioScalar'
        ioScalar particleScalar1[myParticles];

        // Write density rho
        {
            const word what = "rho";
            if (findStrings(r.cloudAttribs_, "rho"))
            {
                Info<< "      dataset " << what << endl;
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    particleScalar1[i++] = pIter().rho();
                }
                writeVariable(varPath/what, &particleScalar1);
            }
        }

        // Write diameter d
        {
            const word what = "d";
            if (findStrings(r.cloudAttribs_, what))
            {
                Info<< "      dataset " << what << endl;
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    particleScalar1[i++] = pIter().d();
                }
                writeVariable(varPath/what, &particleScalar1);
            }
        }

        // Write age
        {
            const word what = "age";
            if (findStrings(r.cloudAttribs_, what))
            {
                Info<< "      dataset " << what << endl;
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    particleScalar1[i++] = pIter().age();
                }
                writeVariable(varPath/what, &particleScalar1);
            }
        }


        // Allocate storage for 3-comp. dataset of type 'ioScalar'
        ioScalar particleScalar3[myParticles*3];

        // Write position
        {
            const word what = "position";
            if (findStrings(r.cloudAttribs_, what))
            {
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    particleScalar3[3*i+0] = pIter().position().x();
                    particleScalar3[3*i+1] = pIter().position().y();
                    particleScalar3[3*i+2] = pIter().position().z();
                    i++;
                }
                writeVariable(varPath/what, &particleScalar3);
            }
        }

        // Write velocity U
        {
            const word what = "U";
            if (findStrings(r.cloudAttribs_, what))
            {
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    particleScalar3[3*i+0] = pIter().U().x();
                    particleScalar3[3*i+1] = pIter().U().y();
                    particleScalar3[3*i+2] = pIter().U().z();
                    i++;
                }
                writeVariable(varPath/what, &particleScalar3);
            }
        }

        // Write slip velocity Us = U - Uc
        {
            const word what = "Us";
            if (findStrings(r.cloudAttribs_, what))
            {
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    particleScalar3[3*i+0] = pIter().U().x() - pIter().Uc().x();
                    particleScalar3[3*i+1] = pIter().U().y() - pIter().Uc().y();
                    particleScalar3[3*i+2] = pIter().U().z() - pIter().Uc().z();
                    i++;
                }
                writeVariable(varPath/what, &particleScalar3);
            }
        }

    }
}


// ************************************************************************* //
