/* cteasFOAM: staged scalar, special-case Couette and M1B-1a laminar
 * pressure-corrected finite-volume modes. The PISO-style projection is
 * initially validated only with a periodic 2D Taylor-Green vortex.
 * No validated arbitrary 3D flow, URANS, HT/HTO, soil or interface yet.
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

    if (momentumMode != "prescribed" && momentumMode != "laminarShear"
        && momentumMode != "pisoLaminar")
    {
        FatalErrorInFunction
            << "Unsupported momentumMode: " << momentumMode
            << ". Valid: prescribed, laminarShear, pisoLaminar"
            << exit(FatalError);
    }

    Foam::dimensionedScalar nu("nu", Foam::dimViscosity, 0.0);
    if (momentumMode == "laminarShear" || momentumMode == "pisoLaminar")
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

    // Only the pressure-corrected mode reads/writes a physical kinematic
    // pressure [m2/s2]. Older scalar and Couette case files need no p field.
    Foam::autoPtr<Foam::volScalarField> pField;
    if (momentumMode == "pisoLaminar")
    {
        pField.reset
        (
            new Foam::volScalarField
            (
                Foam::IOobject
                (
                    "p", runTime.timeName(), mesh,
                    Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }

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

        if (momentumMode == "pisoLaminar")
        {
            // Minimal laminar PISO-style pressure correction. Its full
            // 2-D periodic Taylor-Green reference and flux divergence
            // are independently tested. Full 3-D flow is NOT accepted yet.
            Foam::volScalarField& p = pField();
            Foam::fvVectorMatrix UEqn
            (
                Foam::fvm::ddt(U)
              + Foam::fvm::div(phi, U)
              - Foam::fvm::laplacian(nu, U)
            );
            Foam::solve(UEqn == -Foam::fvc::grad(p));

            Foam::volScalarField rAU("rAU", 1.0/UEqn.A());
            Foam::volVectorField HbyA("HbyA", rAU*UEqn.H());
            Foam::surfaceScalarField phiHbyA
            (
                "phiHbyA",
                Foam::fvc::flux(HbyA)
              + Foam::fvc::interpolate(rAU)*Foam::fvc::ddtCorr(U, phi)
            );
            Foam::fvScalarMatrix pEqn
            (
                Foam::fvm::laplacian(Foam::fvc::interpolate(rAU), p)
             == Foam::fvc::div(phiHbyA)
            );
            // Constant-pressure gauge in the periodic closed domain.
            pEqn.setReference(0, 0.0);
            pEqn.solve();
            phi = phiHbyA - pEqn.flux();
            U = HbyA - rAU*Foam::fvc::grad(p);
            U.correctBoundaryConditions();

            const Foam::volScalarField divPhi(Foam::fvc::div(phi));
            const Foam::scalar maxDivPhi
            (
                Foam::gMax(Foam::mag(divPhi.primitiveField()))
            );
            Foam::scalar netBoundaryFlux = 0.0;
            forAll(phi.boundaryField(), patchI)
            {
                const auto& patchFlux = phi.boundaryField()[patchI];
                forAll(patchFlux, faceI)
                {
                    netBoundaryFlux += patchFlux[faceI];
                }
            }
            Foam::reduce(netBoundaryFlux, Foam::sumOp<Foam::scalar>());
            Foam::Info << "CTEAS_M1B1 max_div_phi=" << maxDivPhi
                       << " boundary_net_flux=" << netBoundaryFlux
                       << " kinetic_energy="
                       << Foam::fvc::domainIntegrate
                          (
                              0.5*Foam::magSqr(U)
                          ).value()
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
               << "; URANS, HT/HTO and soil NOT IMPLEMENTED"
               << Foam::nl;
    return 0;
}
