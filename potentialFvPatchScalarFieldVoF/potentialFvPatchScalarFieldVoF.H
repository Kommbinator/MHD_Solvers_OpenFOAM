#ifndef potentialFvPatchScalarFieldVoF_H
#define potentialFvPatchScalarFieldVoF_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class potentialFvPatchScalarFieldVoF
:
    public mixedFvPatchScalarField
{
    // Private data

public:

    //- Runtime type information
    TypeName("coupledPotentialVoF");


    // Constructors

        //- Construct from patch and internal field
        potentialFvPatchScalarFieldVoF
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        potentialFvPatchScalarFieldVoF
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  potentialFvPatchScalarFieldVoF onto a
        //  new patch
        potentialFvPatchScalarFieldVoF
        (
            const potentialFvPatchScalarFieldVoF&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new potentialFvPatchScalarFieldVoF
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        potentialFvPatchScalarFieldVoF
        (
            const potentialFvPatchScalarFieldVoF&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new potentialFvPatchScalarFieldVoF
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
