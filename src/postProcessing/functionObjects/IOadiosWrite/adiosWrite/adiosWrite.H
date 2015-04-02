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

Class
    Foam::adiosWrite

Description
    Writes selected (vector and scalar) fields and particle clouds to an
    ADIOS file.
    Based on the IOh5Write module created by
       Håkon Strandenes, hakostra@stud.ntnu.no

SourceFiles
    adiosWrite.C
    IOadiosWrite.H

Author
    Norbert Podhorszki, pnorbert@ornl.gov, 2015

\*---------------------------------------------------------------------------*/

#ifndef adiosWrite_H
#define adiosWrite_H

#include "runTimeSelectionTables.H"
#include "interpolation.H"
#include "fvMesh.H"

#include "adios.h"
#include <stdio.h> // sprintf

// If OpenFOAM is compiled with single precision, the IO will be single
// precision.
#if defined(WM_DP)
#   define ADIOS_SCALAR adios_double
    typedef double ioScalar;
#elif defined(WRITE_SP)
#   define ADIOS_SCALAR adios_float
    typedef float ioScalar;
#else
#   error "Neither WM_DP nor WM_SP is defined when compiling adiosWrite.H"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
//class fvMesh;

/*---------------------------------------------------------------------------*\
                         Class adiosWrite Declaration
\*---------------------------------------------------------------------------*/

class adiosWrite
{
protected:   
    // Protected classes

        //- Class template used for grouping field types
        template<class Type>
        class fieldGroup
        :
            public DynamicList<word>
        {
        public:
            //- Construct null
            fieldGroup()
            :
                DynamicList<word>(0)
            {}

        };
    
    
    // Private data

        //- Name of this set of adiosWrite
        word name_;
        
        // Current object registry
        const objectRegistry& obr_;
        
        //- Const reference to fvMesh
        const fvMesh& mesh_;


        // Read from dictionary

            //- Names of volume fields to write
            wordReList objectNames_;
            
            //- Names of clouds to write
            wordReList cloudNames_;
            
            //- Which attributes to write for particle data
            wordReList cloudAttribs_;
            
            //- ADIOS I/O method to use
            word adiosMethod_;

            //- ADIOS I/O method parameters
            string methodParams_;

            //- Write interval (in number of iterations)
            label writeInterval_;
        
        // Number of fields (vector, scalar) to write
        label nFields_;
        
        // Categorized scalar/vector/tensor fields
        fieldGroup<scalar> scalarFields_;
        fieldGroup<vector> vectorFields_;
        
        
        // Used to track when to write
        label nextWrite_;
        label timeSteps_;

        // ADIOS File and other identifiers
        int64_t  fileID_;     // file handler from open()
        int64_t  groupID_;    // group of output variables
        uint64_t outputSize_; // number of bytes we write in one timestep by the process
        
        // Array containing number of cells per process (often used list)
        List<label> nCells_;

        // number of data points in output Cell array on each process
        //   calculated in meshDefineCells, used in meshWriteCells
        //   only the local size is used currently
        List<label> cellDataSizes_;
            
        // MPI parallel code variables
        MPI_Comm comm_; // Pstream:: does not provide communicator
        // rank = Pstream::myProcNo()
        // number of processors = Pstream::nProcs()
        // process id: Pstream::procID (int procNo), 1 for rank 0
            

    // Private Member Functions
      
      // General functions
      
          // Create and open dataset for all data
          void open();
          
          // Close dataset
          void close();
          

      
      // Functions for handling field data
          
          // Field define wrapper
          void fieldDefine();

          // Write scalar data
          void fieldDefineScalar();

          // Write vector data
          void fieldDefineVector();

          // Field write wrapper
          void fieldWrite();
          
          // Write scalar data
          void fieldWriteScalar();
          
          // Write vector data
          void fieldWriteVector();
      
      
      // Functions for handling particle cloud data
          
          // Cloud define wrapper
          void cloudDefine();

          // Cloud write wrapper
          void cloudWrite();
          
          // Write attributes for cloud particles
          void cloudWriteAttrib
          (
              label myParticles,
              label offset,
              label nTot,
              label nCmps,
              void* databuf,
              char* datasetName,
              ioScalar atype
          );
          
      // Functions for handling mesh 
      
          // Mesh define wrapper
          void meshDefine();

          // Define points
          void meshDefinePoints();

          // Define cells
          void meshDefineCells();

          // Define Boundary meshes
          void meshDefineBoundaries();

          // Mesh write wrapper
          void meshWrite();
          
          // Write points
          void meshWritePoints();
          
          // Write cells
          void meshWriteCells();
          
      // Miscallenous definitions

          //- Append fieldName to the appropriate group
          label appendFieldGroup(const word& fieldName, const word& fieldType);

          //- Classify field types, returns the number of fields
          label classifyFields();
          
          //- Disallow default bitwise copy construct
          adiosWrite(const adiosWrite&);

          //- Disallow default bitwise assignment
          void operator=(const adiosWrite&);
      
          // define ADIOS variables (before writing them)
          void defineVars ();
          void deleteDefinitions();
      
      //- Pointers to cell shape models
      static const cellModel* unknownModel;
      static const cellModel* tetModel;
      static const cellModel* pyrModel;
      static const cellModel* prismModel;
      static const cellModel* hexModel;


public:

    //- Runtime type information
    TypeName("adiosWrite");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        adiosWrite
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~adiosWrite();


    // Member Functions

        //- Return name of the adiosWrite
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the adiosWrite data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();
        
        //- Called when time was set at the end of the Time::operator++,
        //- currently does nothing
        virtual void timeSet();

        //- Write a timestep to file
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
