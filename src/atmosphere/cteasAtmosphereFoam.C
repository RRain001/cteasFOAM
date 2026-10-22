/* cteasFOAM M1A: prescribed, divergence-free U; conservative passive-scalar FVM.
 * No momentum/URANS closure, HT/HTO conversion, deposition, or soil coupling.
 * The scalar c is dimensionless; activity units are deliberately deferred.
 */
#include "fvCFD.H"

int main(int argc, char* argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Foam::IOdictionary transportProperties
    (
        Foam::IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE
        )
    );

    const Foam::dimensionedScalar D
    (
        "D",
        Foam::dimViscosity,
        transportProperties
    );

    if (D.value() < 0)
    {
        FatalErrorInFunction
            << "Passive-scalar diffusivity D must be non-negative"
            << exit(FatalError);
    }

    Foam::volVectorField U
    (
        Foam::IOobject
        (
            "U", runTime.timeName(), mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE
        ),
        mesh
    );

    Foam::volScalarField c
    (
        Foam::IOobject
        (
            "c", runTime.timeName(), mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Incompressible volumetric face flux [m3/s] from a prescribed velocity.
    // Valid conservation tests require div(phi)=0, which is the case
    // for the closed, zero-velocity M1A diffusion benchmark.
    Foam::surfaceScalarField phi("phi", Foam::fvc::flux(U));

    Foam::Info << "M1A passive scalar; prescribed velocity; D="
               << D.value() << " m2/s" << Foam::nl
               << "CTEAS_M1A initial_integral="
               << Foam::fvc::domainIntegrate(c).value() << Foam::nl;

    while (runTime.loop())
    {
        Foam::Info << "Time = " << runTime.timeName() << Foam::nl;
        Foam::fvScalarMatrix cEqn
        (
            Foam::fvm::ddt(c)
          + Foam::fvm::div(phi, c)
          - Foam::fvm::laplacian(D, c)
        );
        cEqn.solve();

        Foam::Info << "CTEAS_M1A integral="
                   << Foam::fvc::domainIntegrate(c).value() << Foam::nl;
        runTime.write();
    }

    Foam::Info << "M1A completed; atmospheric URANS/HT/HTO NOT IMPLEMENTED"
               << Foam::nl;
    return 0;
}
