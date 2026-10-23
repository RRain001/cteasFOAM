/* cteasFOAM staged atmospheric solver: M1A scalar + M1B-0 laminar shear.
 * M1B-0 solves viscous transient momentum diffusion ONLY in an exactly
 * divergence-free planar shear flow. It is not a general pressure-projection,
 * Navier-Stokes, buoyant, turbulent or URANS model.
 * HT/HTO activity, deposition and air/soil coupling remain unimplemented.
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
            "transportProperties", runTime.constant(), mesh,
            Foam::IOobject::MUST_READ, Foam::IOobject::NO_WRITE
        )
    );

    const Foam::dimensionedScalar D
    (
        "D", Foam::dimViscosity, transportProperties
    );

    if (D.value() < 0)
    {
        FatalErrorInFunction
            << "Passive scalar diffusivity D must be non-negative"
            << exit(FatalError);
    }

    // Preserve original M1A behavior unless an explicit experimental mode
    // is selected. Do not imply that a solved planar shear velocity is a
    // general incompressible pressure-velocity coupled wind field.
    const Foam::word momentumMode
    (
        transportProperties.lookupOrDefault<Foam::word>
        (
            "momentumMode", "prescribed"
        )
    );

    if (momentumMode != "prescribed" && momentumMode != "laminarShear")
    {
        FatalErrorInFunction
            << "Unsupported momentumMode: " << momentumMode
            << ". Valid: prescribed, laminarShear"
            << exit(FatalError);
    }

    Foam::dimensionedScalar nu("nu", Foam::dimViscosity, 0.0);
    if (momentumMode == "laminarShear")
    {
        nu = Foam::dimensionedScalar
        (
            "nu", Foam::dimViscosity, transportProperties
        );
        if (nu.value() <= 0)
        {
            FatalErrorInFunction
                << "Laminar momentum viscosity nu must be strictly positive"
                << exit(FatalError);
        }
    }

    Foam::volVectorField U
    (
        Foam::IOobject
        (
            "U", runTime.timeName(), mesh,
            Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE
        ),
        mesh
    );

    Foam::volScalarField c
    (
        Foam::IOobject
        (
            "c", runTime.timeName(), mesh,
            Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE
        ),
        mesh
    );

    Foam::surfaceScalarField phi("phi", Foam::fvc::flux(U));
    Foam::Info << "cteasFOAM mode=" << momentumMode
               << " D=" << D.value() << Foam::nl
               << "CTEAS_M1A initial_integral="
               << Foam::fvc::domainIntegrate(c).value() << Foam::nl;

    while (runTime.loop())
    {
        Foam::Info << "Time = " << runTime.timeName() << Foam::nl;

        if (momentumMode == "laminarShear")
        {
            // For planar Couette: U=(u(y,t),0,0). Its divergence is identically
            // zero, both pressure gradient and nonlinear convection vanish.
            // This is a restricted momentum PDE, NOT a 3-D URANS solver.
            Foam::fvVectorMatrix UEqn
            (
                Foam::fvm::ddt(U)
              - Foam::fvm::laplacian(nu, U)
            );
            UEqn.solve();
            phi = Foam::fvc::flux(U);
            Foam::Info << "CTEAS_M1B momentum_integral="
                       << Foam::fvc::domainIntegrate(U).value()
                       << Foam::nl;
        }

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

    Foam::Info << "Completed mode=" << momentumMode
               << "; full URANS, HT/HTO and soil NOT IMPLEMENTED"
               << Foam::nl;
    return 0;
}
