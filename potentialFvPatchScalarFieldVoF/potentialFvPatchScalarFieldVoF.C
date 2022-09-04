#include "potentialFvPatchScalarFieldVoF.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

potentialFvPatchScalarFieldVoF::
potentialFvPatchScalarFieldVoF
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


potentialFvPatchScalarFieldVoF::
potentialFvPatchScalarFieldVoF
(
    const potentialFvPatchScalarFieldVoF& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


potentialFvPatchScalarFieldVoF::
potentialFvPatchScalarFieldVoF
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "potentialFvPatchScalarFieldVoF::"
            "potentialFvPatchScalarFieldVoF\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


potentialFvPatchScalarFieldVoF::
potentialFvPatchScalarFieldVoF
(
    const potentialFvPatchScalarFieldVoF& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void potentialFvPatchScalarFieldVoF::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());


    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchI = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchI];

    tmp<scalarField> intFld = patchInternalField();

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const fvPatchScalarField& nbrField =  // auf dem Rand
    refCast<const fvPatchScalarField>
    (
        nbrPatch.lookupPatchField<volScalarField, scalar>
        (
             this->internalField().name()
        )
    );
    // Swap to obtain full local values of neighbour internal field
    scalarField nbrIntFld(nbrField.patchInternalField());
    mpp.distribute(nbrIntFld); // verteilt das korrekt zwischen den Prozessoren, wichtig!!!

    const fvPatchField<scalar>& sigmaEMfield = patch().lookupPatchField<volScalarField, scalar>("sigmaEMfield");
    const fvPatchField<scalar>& nbrSigma = nbrField.patch().lookupPatchField<volScalarField, scalar>("sigmaEMfield");

    // Swap to obtain full local values of neighbour kappa*delta
    scalarField nbrSigmaDelta(nbrSigma*nbrPatch.deltaCoeffs());
    mpp.distribute(nbrSigmaDelta); // das ist parallel sehr wichtig!

    tmp<scalarField> sigmaDelta = sigmaEMfield*patch().deltaCoeffs();


// *********************************************************************
// *********************************************************************
// *********************************************************************

    this->refValue() = nbrIntFld;
    this->refGrad() = 0.0;

    valueFraction() = nbrSigmaDelta /
	              (nbrSigmaDelta + sigmaDelta);


// *********************************************************************
// *********************************************************************
// *********************************************************************

    mixedFvPatchScalarField::updateCoeffs();
    // Restore tag
    UPstream::msgType() = oldTag;
}


void potentialFvPatchScalarFieldVoF::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    potentialFvPatchScalarFieldVoF
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
