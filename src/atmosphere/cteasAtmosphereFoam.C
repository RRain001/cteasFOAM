/* cteasFOAM: staged scalar, special-case Couette and M1B-1a laminar
 * pressure-corrected finite-volume modes. The PISO-style projection is
 * initially validated only with a periodic 2D Taylor-Green vortex.
 * No validated arbitrary 3D flow, URANS, HT/HTO, soil or interface yet.
 */
#include "fvCFD.H"
#include "constrainHbyA.H"
#include "adjustPhi.H"
#include "constrainPressure.H"
#include "processorCyclicPolyPatch.H"
#include <cmath>
#include <sstream>
#include <iomanip>

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

    // Optional M1C-0 two-species activity prototype. Default none keeps all
    // existing dimensionless passive-scalar cases exactly unchanged.
    const Foam::word activityMode
    (
        transportProperties.lookupOrDefault<Foam::word>
        (
            "activityMode", "none"
        )
    );
    if (activityMode != "none" && activityMode != "htHto")
    {
        FatalErrorInFunction
            << "Unsupported activityMode: " << activityMode
            << ". Valid: none, htHto" << exit(FatalError);
    }

    // M1C-17 is a fail-closed diagnostic, not a positivity correction.
    // Disabled by default to leave all frozen legacy equations unchanged.
    const Foam::Switch enforceNonnegativeActivity
    (
        transportProperties.lookupOrDefault<Foam::Switch>
        (
            "enforceNonnegativeActivity", false
        )
    );
    if (enforceNonnegativeActivity && activityMode != "htHto")
    {
        FatalErrorInFunction
            << "enforceNonnegativeActivity requires activityMode htHto"
            << exit(FatalError);
    }

    // M1C-18 restricted high/low finite-volume convex blending.
    // Not face-local FCT or MULES; default preserves the existing solver.
    const Foam::word activityConvectionBlend
    (
        transportProperties.lookupOrDefault<Foam::word>
        (
            "activityConvectionBlend", "none"
        )
    );
    if
    (
        activityConvectionBlend != "none"
     && activityConvectionBlend != "globalPositive"
    )
    {
        FatalErrorInFunction
            << "Unsupported activityConvectionBlend: "
            << activityConvectionBlend << exit(FatalError);
    }
    if (activityConvectionBlend == "globalPositive" && activityMode != "htHto")
    {
        FatalErrorInFunction
            << "globalPositive requires activityMode htHto"
            << exit(FatalError);
    }

    Foam::dimensionedScalar kOx
    (
        "kOx", Foam::dimless/Foam::dimTime, 0.0
    );
    Foam::dimensionedScalar lambdaT
    (
        "lambdaT", Foam::dimless/Foam::dimTime, 0.0
    );
    Foam::dimensionedScalar DHT("DHT", Foam::dimViscosity, 0.0);
    Foam::dimensionedScalar DHTO("DHTO", Foam::dimViscosity, 0.0);
    if (activityMode == "htHto")
    {
        kOx = Foam::dimensionedScalar
        (
            "kOx", Foam::dimless/Foam::dimTime, transportProperties
        );
        lambdaT = Foam::dimensionedScalar
        (
            "lambdaT", Foam::dimless/Foam::dimTime, transportProperties
        );
        DHT = Foam::dimensionedScalar
        (
            "DHT", Foam::dimViscosity, transportProperties
        );
        DHTO = Foam::dimensionedScalar
        (
            "DHTO", Foam::dimViscosity, transportProperties
        );
        if (kOx.value() < 0 || lambdaT.value() < 0)
        {
            FatalErrorInFunction
                << "M1C reaction rates kOx and lambdaT must be non-negative"
                << exit(FatalError);
        }
        if (DHT.value() < 0 || DHTO.value() < 0)
        {
            FatalErrorInFunction
                << "M1C diffusivities DHT and DHTO must be non-negative"
                << exit(FatalError);
        }
    }

    // Optional M1C-5 finite-time Bq/s injection and M1C-6 normalized
    // 3-D spatial kernel. Emission rates are integrated Bq/s (s^-2), not
    // spatial source densities. Neither mode is a physical stack model.
    const Foam::dimensionSet emissionDims(0, 0, -2, 0, 0, 0, 0);
    Foam::dimensionedScalar sourceHTRate
    (
        "sourceHT_Bq_s", emissionDims, 0.0
    );
    Foam::dimensionedScalar sourceHTORate
    (
        "sourceHTO_Bq_s", emissionDims, 0.0
    );
    Foam::dimensionedScalar sourceStart("sourceStart", Foam::dimTime, 0.0);
    Foam::dimensionedScalar sourceEnd("sourceEnd", Foam::dimTime, 0.0);
    const bool hasActivitySource =
        activityMode == "htHto"
     && (transportProperties.found("sourceHT_Bq_s")
      || transportProperties.found("sourceHTO_Bq_s"));
    const Foam::word sourceShape
    (
        transportProperties.lookupOrDefault<Foam::word>
        (
            "sourceShape", "uniform"
        )
    );
    if
    (
        transportProperties.found("sourceShape") && !hasActivitySource
    )
    {
        FatalErrorInFunction
            << "M1C sourceShape requires configured Bq/s emission rates"
            << exit(FatalError);
    }
    if (sourceShape != "uniform" && sourceShape != "gaussian")
    {
        FatalErrorInFunction
            << "M1C unsupported sourceShape: " << sourceShape
            << ". Valid: uniform, gaussian" << exit(FatalError);
    }
    // M1C-15: opt-in conservative polyhedral Gaussian volume quadrature.
    // The original cell-centre kernel remains the default, exactly.
    const Foam::word sourceSpatialQuadrature
    (
        transportProperties.lookupOrDefault<Foam::word>
        (
            "sourceSpatialQuadrature", "cellCentre"
        )
    );
    if
    (
        sourceSpatialQuadrature != "cellCentre"
     && sourceSpatialQuadrature != "cellTet4"
    )
    {
        FatalErrorInFunction
            << "M1C unsupported sourceSpatialQuadrature: "
            << sourceSpatialQuadrature << ". Valid: cellCentre, cellTet4"
            << exit(FatalError);
    }
    const Foam::Switch writeSourceProfile
    (
        transportProperties.lookupOrDefault<Foam::Switch>
        (
            "writeSourceProfile", false
        )
    );
    if
    (
        (transportProperties.found("sourceSpatialQuadrature")
      || transportProperties.found("writeSourceProfile"))
     && (!hasActivitySource || sourceShape != "gaussian")
    )
    {
        FatalErrorInFunction
            << "M1C sourceSpatialQuadrature and writeSourceProfile "
            << "require a configured Gaussian activity source"
            << exit(FatalError);
    }

    // Dimensionless spatial profile whose volume-weighted mean is 1.
    // Multiplying the existing rate/domain-volume density by this field
    // conserves exactly the requested Bq/s across native MPI cells.
    Foam::volScalarField sourceProfile
    (
        Foam::IOobject
        (
            "sourceProfile", runTime.timeName(), mesh,
            Foam::IOobject::NO_READ, Foam::IOobject::NO_WRITE
        ),
        mesh,
        Foam::dimensionedScalar("one", Foam::dimless, 1.0)
    );
    Foam::scalar sourceDomainVolume = 1.0;
    if (hasActivitySource)
    {
        if
        (
            !transportProperties.found("sourceStart")
         || !transportProperties.found("sourceEnd")
        )
        {
            FatalErrorInFunction
                << "M1C activity emission requires sourceStart and sourceEnd"
                << exit(FatalError);
        }
        if (transportProperties.found("sourceHT_Bq_s"))
        {
            sourceHTRate = Foam::dimensionedScalar
            (
                "sourceHT_Bq_s", emissionDims, transportProperties
            );
        }
        if (transportProperties.found("sourceHTO_Bq_s"))
        {
            sourceHTORate = Foam::dimensionedScalar
            (
                "sourceHTO_Bq_s", emissionDims, transportProperties
            );
        }
        sourceStart = Foam::dimensionedScalar
        (
            "sourceStart", Foam::dimTime, transportProperties
        );
        sourceEnd = Foam::dimensionedScalar
        (
            "sourceEnd", Foam::dimTime, transportProperties
        );
        if
        (
            !std::isfinite(sourceHTRate.value())
         || !std::isfinite(sourceHTORate.value())
         || !std::isfinite(sourceStart.value())
         || !std::isfinite(sourceEnd.value())
         || sourceHTRate.value() < 0 || sourceHTORate.value() < 0
         || sourceStart.value() < 0
         || sourceEnd.value() <= sourceStart.value()
        )
        {
            FatalErrorInFunction
                << "M1C activity emission requires nonnegative finite Bq/s "
                << "rates and finite sourceEnd > sourceStart >= 0"
                << exit(FatalError);
        }
        // gSum is global across all MPI subdomains, so the summed emission
        // over the physical domain is exactly the prescribed total Bq/s.
        sourceDomainVolume = Foam::gSum(mesh.V());
        if (!std::isfinite(sourceDomainVolume) || sourceDomainVolume <= 0)
        {
            FatalErrorInFunction
                << "M1C activity emission requires positive finite volume"
                << exit(FatalError);
        }
        if (sourceShape == "gaussian")
        {
            if
            (
                !transportProperties.found("sourceCenter")
             || !transportProperties.found("sourceSigma")
            )
            {
                FatalErrorInFunction
                    << "M1C gaussian source requires sourceCenter and sourceSigma"
                    << exit(FatalError);
            }
            const Foam::dimensionedVector sourceCenter
            (
                "sourceCenter", Foam::dimLength, transportProperties
            );
            const Foam::dimensionedScalar sourceSigma
            (
                "sourceSigma", Foam::dimLength, transportProperties
            );
            const Foam::vector center = sourceCenter.value();
            const Foam::scalar sigma = sourceSigma.value();
            if
            (
                !std::isfinite(center.x())
             || !std::isfinite(center.y())
             || !std::isfinite(center.z())
             || !std::isfinite(sigma) || sigma <= 0
            )
            {
                FatalErrorInFunction
                    << "M1C gaussian source requires finite center and positive sigma"
                    << exit(FatalError);
            }
            const Foam::scalar invTwoSigma2 = 1.0/(2*sigma*sigma);
            const Foam::scalar tetA =
                (5.0 + 3.0*std::sqrt(5.0))/20.0;
            const Foam::scalar tetB =
                (5.0 - std::sqrt(5.0))/20.0;
            const Foam::faceList& faces = mesh.faces();
            const Foam::pointField& points = mesh.points();
            const Foam::vectorField& faceCenters = mesh.faceCentres();
            const Foam::cellList& cells = mesh.cells();
            forAll(sourceProfile.primitiveField(), cellI)
            {
                if (sourceSpatialQuadrature == "cellCentre")
                {
                    sourceProfile.primitiveFieldRef()[cellI] = std::exp
                    (
                        -Foam::magSqr(mesh.C()[cellI] - center)
                        *invTwoSigma2
                    );
                    continue;
                }

                // Split each true native polyhedral face about its centre
                // into triangles and form positive-volume tetrahedra
                // with the actual native cell volume centre.
                const Foam::vector apex = mesh.C()[cellI];
                Foam::scalar integral = 0.0;
                Foam::scalar tetraVolumeSum = 0.0;
                const Foam::cell& cellFaces = cells[cellI];
                forAll(cellFaces, faceLocalI)
                {
                    const Foam::label faceI = cellFaces[faceLocalI];
                    const Foam::face& f = faces[faceI];
                    if (f.size() < 3)
                    {
                        FatalErrorInFunction
                            << "M1C cellTet4 requires valid polygonal faces"
                            << exit(FatalError);
                    }
                    const Foam::vector faceCenter = faceCenters[faceI];
                    forAll(f, vertexI)
                    {
                        const Foam::vector triangleA = points[f[vertexI]];
                        const Foam::vector triangleB =
                            points[f[(vertexI + 1) % f.size()]];
                        const Foam::scalar tetraVolume =
                            Foam::mag
                            (
                                (faceCenter - apex)
                              & ((triangleA - apex)
                              ^ (triangleB - apex))
                            )/6.0;
                        if
                        (
                            !std::isfinite(tetraVolume)
                         || tetraVolume <= 0.0
                        )
                        {
                            FatalErrorInFunction
                                << "M1C cellTet4 invalid native tetra volume"
                                << exit(FatalError);
                        }
                        tetraVolumeSum += tetraVolume;
                        const Foam::vector tet[4] =
                        {
                            apex, faceCenter, triangleA, triangleB
                        };
                        const Foam::vector cornerSum =
                            tet[0] + tet[1] + tet[2] + tet[3];
                        for (Foam::label q = 0; q < 4; ++q)
                        {
                            const Foam::vector sample =
                                tetB*cornerSum + (tetA-tetB)*tet[q];
                            const Foam::scalar kernel = std::exp
                            (
                                -Foam::magSqr(sample-center)*invTwoSigma2
                            );
                            if (!std::isfinite(kernel))
                            {
                                FatalErrorInFunction
                                    << "M1C cellTet4 nonfinite kernel"
                                    << exit(FatalError);
                            }
                            integral += tetraVolume*kernel/4.0;
                        }
                    }
                }
                if
                (
                    !std::isfinite(integral) || integral <= 0.0
                 || !std::isfinite(tetraVolumeSum)
                 || Foam::mag(tetraVolumeSum/mesh.V()[cellI]-1.0)
                       > 0.005
                )
                {
                    FatalErrorInFunction
                        << "M1C cellTet4 invalid native volume quadrature"
                        << " cell=" << cellI
                        << " tetraVolume=" << tetraVolumeSum
                        << " realVolume=" << mesh.V()[cellI]
                        << " gaussianIntegral=" << integral
                        << exit(FatalError);
                }
                sourceProfile.primitiveFieldRef()[cellI] =
                    integral/mesh.V()[cellI];
            }
            const Foam::scalar weightedVolume = Foam::gSum
            (
                mesh.V()*sourceProfile.primitiveField()
            );
            if (!std::isfinite(weightedVolume) || weightedVolume <= 0)
            {
                FatalErrorInFunction
                    << "M1C gaussian source has zero/nonfinite global support"
                    << exit(FatalError);
            }
            sourceProfile.primitiveFieldRef() *=
                sourceDomainVolume/weightedVolume;
            Foam::Info << "CTEAS_M1C sourceShape=gaussian"
                       << " weighted_volume_m3=" << weightedVolume
                       << " domain_volume_m3=" << sourceDomainVolume
                       << " sourceSpatialQuadrature="
                       << sourceSpatialQuadrature
                       << Foam::nl;
            if (writeSourceProfile)
            {
                if (!sourceProfile.write())
                {
                    FatalErrorInFunction
                        << "M1C requested sourceProfile diagnostic write failed"
                        << exit(FatalError);
                }
            }
        }
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

    // A constant acceleration has the same streamwise momentum effect as
    // a prescribed uniform pressure gradient in a periodic channel.
    // Default zero preserves every earlier M1A/M1B benchmark unchanged.
    Foam::dimensionedVector drivingAcceleration
    (
        "drivingAcceleration", Foam::dimAcceleration, Foam::vector::zero
    );
    if
    (
        momentumMode == "pisoLaminar"
     && transportProperties.found("drivingAcceleration")
    )
    {
        drivingAcceleration = Foam::dimensionedVector
        (
            "drivingAcceleration", Foam::dimAcceleration,
            transportProperties
        );
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

    Foam::autoPtr<Foam::volScalarField> htField;
    Foam::autoPtr<Foam::volScalarField> htoField;
    if (activityMode == "htHto")
    {
        htField.reset
        (
            new Foam::volScalarField
            (
                Foam::IOobject
                (
                    "HT", runTime.timeName(), mesh,
                    Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
        htoField.reset
        (
            new Foam::volScalarField
            (
                Foam::IOobject
                (
                    "HTO", runTime.timeName(), mesh,
                    Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
        const Foam::dimensionSet activityDims(0, -3, -1, 0, 0, 0, 0);
        if
        (
            htField().dimensions() != activityDims
         || htoField().dimensions() != activityDims
        )
        {
            FatalErrorInFunction
                << "M1C HT/HTO fields require Bq/m3 dimensions "
                << activityDims << exit(FatalError);
        }
    }

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

    if (activityMode == "htHto")
    {
        Foam::Info.stream().precision(17);
    }
    Foam::surfaceScalarField phi("phi", Foam::fvc::flux(U));
    Foam::Info << "cteasFOAM mode=" << momentumMode
               << " D=" << D.value() << Foam::nl
               << "CTEAS_M1A initial_integral="
               << Foam::fvc::domainIntegrate(c).value() << Foam::nl;

    if (activityConvectionBlend == "globalPositive")
    {
        // This first research baseline is serial, prescribed-U, cyclic-only.
        // In particular, do not infer MPI-coupled flux conservation here.
        if (Foam::Pstream::parRun() || momentumMode != "prescribed")
        {
            FatalErrorInFunction
                << "M1C18_GLOBAL_BLEND_UNSUPPORTED requires serial prescribed flow"
                << exit(FatalError);
        }
        forAll(mesh.boundaryMesh(), patchI)
        {
            if (mesh.boundaryMesh()[patchI].type() != "cyclic")
            {
                FatalErrorInFunction
                    << "M1C18_GLOBAL_BLEND_UNSUPPORTED requires cyclic-only boundary; patch="
                    << mesh.boundaryMesh()[patchI].name()
                    << exit(FatalError);
            }
        }
    }

    // Use global extrema on all MPI ranks, including at restart/time zero.
    // The tolerance is detection-only: NO clipping or activity correction.
    const auto requireNonnegativeActivity = [&](const char* phase)
    {
        if (!enforceNonnegativeActivity) return;
        // Scalar min/max reductions alone can silently ignore NaN on some
        // implementations. Count nonfinite cells explicitly, then reduce
        // before any rank may abort (including on a decomposed mesh).
        Foam::label htNonfinite = 0, htoNonfinite = 0;
        const Foam::scalarField& htCells = htField().primitiveField();
        const Foam::scalarField& htoCells = htoField().primitiveField();
        forAll(htCells, cellI)
        {
            if (!std::isfinite(htCells[cellI])) ++htNonfinite;
        }
        forAll(htoCells, cellI)
        {
            if (!std::isfinite(htoCells[cellI])) ++htoNonfinite;
        }
        Foam::reduce(htNonfinite, Foam::sumOp<Foam::label>());
        Foam::reduce(htoNonfinite, Foam::sumOp<Foam::label>());
        const Foam::scalar htMin = Foam::gMin(htField().primitiveField());
        const Foam::scalar htMax = Foam::gMax(htField().primitiveField());
        const Foam::scalar htoMin = Foam::gMin(htoField().primitiveField());
        const Foam::scalar htoMax = Foam::gMax(htoField().primitiveField());
        const Foam::scalar tol = -1e-8;
        if (htNonfinite || !std::isfinite(htMin) || !std::isfinite(htMax) || htMin < tol)
        {
            FatalErrorInFunction
                << "CTEAS_M1C_POSITIVITY_FAILURE species=HT phase=" << phase
                << " minimum_Bq_m3=" << htMin << " threshold_Bq_m3=" << tol
                << " nonfinite_cell_count=" << htNonfinite
                << exit(FatalError);
        }
        if (htoNonfinite || !std::isfinite(htoMin) || !std::isfinite(htoMax) || htoMin < tol)
        {
            FatalErrorInFunction
                << "CTEAS_M1C_POSITIVITY_FAILURE species=HTO phase=" << phase
                << " minimum_Bq_m3=" << htoMin << " threshold_Bq_m3=" << tol
                << " nonfinite_cell_count=" << htoNonfinite
                << exit(FatalError);
        }
    };
    if (activityMode == "htHto") requireNonnegativeActivity("initial");

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
            Foam::solve(UEqn == -Foam::fvc::grad(p) + drivingAcceleration);

            // Multiple PISO pressure correctors are optional and apply
            // only to this pressure-corrected laminar mode. Default 1
            // retains all previously validated Taylor-Green behaviour.
            const Foam::label nPressureCorrectors =
                transportProperties.lookupOrDefault<Foam::label>
                (
                    "nPressureCorrectors", 1
                );
            if (nPressureCorrectors < 1 || nPressureCorrectors > 12)
            {
                FatalErrorInFunction
                    << "nPressureCorrectors must lie in [1, 12]"
                    << exit(FatalError);
            }
            // Inner non-orthogonal correction is distinct from the
            // existing outer PISO pressure-velocity corrector count.
            // Zero (one p solve) preserves all validated older cases.
            const Foam::label nNonOrthogonalCorrectors =
                transportProperties.lookupOrDefault<Foam::label>
                (
                    "nNonOrthogonalCorrectors", 0
                );
            if
            (
                nNonOrthogonalCorrectors < 0
             || nNonOrthogonalCorrectors > 6
            )
            {
                FatalErrorInFunction
                    << "nNonOrthogonalCorrectors must lie in [0, 6]"
                    << exit(FatalError);
            }
            Foam::Info << "CTEAS_M1B1 nonorthogonal_correctors="
                       << nNonOrthogonalCorrectors << Foam::nl;
            Foam::volScalarField rAU("rAU", 1.0/UEqn.A());
            for
            (
                Foam::label pressureCorrector = 0;
                pressureCorrector < nPressureCorrectors;
                ++pressureCorrector
            )
            {
                // Retain explicit body-force acceleration in the corrected
                // momentum source as well as in the predictor equation.
                // Without this term PISO projection silently drops the drive.
                Foam::volVectorField HbyA
                (
                    "HbyA",
                    Foam::constrainHbyA
                    (
                        rAU*(UEqn.H() + drivingAcceleration), U, p
                    )
                );
                Foam::surfaceScalarField phiHbyA
                (
                    "phiHbyA",
                    Foam::fvc::flux(HbyA)
                  + Foam::fvc::interpolate(rAU)*Foam::fvc::ddtCorr(U, phi)
                );
                // Standard pressure/flux boundary consistency for fixed-U
                // no-slip walls; invisible in the original fully cyclic TGV.
                Foam::adjustPhi(phiHbyA, U, p);
                Foam::constrainPressure(p, U, phiHbyA, rAU);
                // Rebuild the pressure equation for each non-orthogonal
                // correction because corrected snGrad/laplacian operators
                // contain explicit pressure contributions. Correct phi
                // ONLY on the final inner pass (OpenFOAM PISO convention).
                // All ranks enter setReference for each pressure solve.
                for
                (
                    Foam::label nonOrthI = 0;
                    nonOrthI <= nNonOrthogonalCorrectors;
                    ++nonOrthI
                )
                {
                    Foam::fvScalarMatrix pEqn
                    (
                        Foam::fvm::laplacian
                        (
                            Foam::fvc::interpolate(rAU), p
                        )
                     == Foam::fvc::div(phiHbyA)
                    );
                    const Foam::label pRefCell =
                        Foam::Pstream::master() ? 0 : -1;
                    pEqn.setReference(pRefCell, 0.0);
                    pEqn.solve();
                    if (nonOrthI == nNonOrthogonalCorrectors)
                    {
                        phi = phiHbyA - pEqn.flux();
                    }
                }
                U = HbyA - rAU*Foam::fvc::grad(p);
                U.correctBoundaryConditions();
            }

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
            // Optional direct measurement of true open-boundary flow.
            // A net zero flux alone cannot prove nonzero inlet/outlet
            // throughflow, so audit each named pressure patch separately.
            if
            (
                transportProperties.lookupOrDefault<Foam::Switch>
                (
                    "auditOpenPressureFlux", false
                )
            )
            {
                Foam::scalar qIn = 0.0;
                Foam::scalar qOut = 0.0;
                const Foam::label inPatch =
                    mesh.boundaryMesh().findPatchID("xMin");
                const Foam::label outPatch =
                    mesh.boundaryMesh().findPatchID("xMax");
                if (inPatch >= 0)
                {
                    qIn = Foam::sum(phi.boundaryField()[inPatch]);
                }
                if (outPatch >= 0)
                {
                    qOut = Foam::sum(phi.boundaryField()[outPatch]);
                }
                Foam::reduce(qIn, Foam::sumOp<Foam::scalar>());
                Foam::reduce(qOut, Foam::sumOp<Foam::scalar>());
                Foam::Info << "CTEAS_M1B1 pressure_inlet_flux=" << qIn
                           << " pressure_outlet_flux=" << qOut
                           << Foam::nl;
            }
            Foam::Info << "CTEAS_M1B1 max_div_phi=" << maxDivPhi
                       << " boundary_net_flux=" << netBoundaryFlux
                       << " kinetic_energy="
                       << Foam::fvc::domainIntegrate
                          (
                              0.5*Foam::magSqr(U)
                          ).value()
                       << Foam::nl;
        }

        // Optional equation-level audit for real nonuniform scalar inflow.
        // Unlike the historical uniform c=1 closed-volume regression, this
        // uses the *same implicit FVM operator* used to advance c, including
        // its actual selected upwind boundary coefficients. c is
        // dimensionless; its integrated content is measured in m3, NOT Bq.
        const Foam::Switch auditPassiveScalarFlux
        (
            transportProperties.lookupOrDefault<Foam::Switch>
            (
                "auditPassiveScalarFlux", false
            )
        );
        const Foam::scalar oldScalarInventory =
            auditPassiveScalarFlux
          ? Foam::fvc::domainIntegrate(c).value()
          : 0.0;
        Foam::fvScalarMatrix cEqn
        (
            Foam::fvm::ddt(c)
          + Foam::fvm::div(phi, c)
          - Foam::fvm::laplacian(D, c)
        );
        cEqn.solve();

        if (auditPassiveScalarFlux)
        {
            // fluxRequired(c) is explicitly enabled for this experiment.
            // Matrix flux is oriented outward at boundary patches and
            // includes the flux discretization actually solved by cEqn.
            const Foam::surfaceScalarField cFaceFlux(cEqn.flux());
            Foam::scalar cFluxNet = 0.0;
            Foam::scalar cFluxXMin = 0.0;
            Foam::scalar cFluxXMax = 0.0;
            const Foam::label xMinPatch =
                mesh.boundaryMesh().findPatchID("xMin");
            const Foam::label xMaxPatch =
                mesh.boundaryMesh().findPatchID("xMax");
            if (xMinPatch < 0 || xMaxPatch < 0)
            {
                FatalErrorInFunction
                    << "auditPassiveScalarFlux needs xMin/xMax patches"
                    << exit(FatalError);
            }
            forAll(cFaceFlux.boundaryField(), patchI)
            {
                const Foam::scalar patchSum =
                    Foam::sum(cFaceFlux.boundaryField()[patchI]);
                cFluxNet += patchSum;
                if (patchI == xMinPatch) cFluxXMin = patchSum;
                if (patchI == xMaxPatch) cFluxXMax = patchSum;
            }
            Foam::reduce(cFluxNet, Foam::sumOp<Foam::scalar>());
            Foam::reduce(cFluxXMin, Foam::sumOp<Foam::scalar>());
            Foam::reduce(cFluxXMax, Foam::sumOp<Foam::scalar>());
            const Foam::scalar newScalarInventory =
                Foam::fvc::domainIntegrate(c).value();
            const Foam::scalar cBudgetResidual =
                newScalarInventory - oldScalarInventory
              + runTime.deltaTValue()*cFluxNet;
            Foam::Info << "CTEAS_M1A scalar_flux_xMin=" << cFluxXMin
                       << " scalar_flux_xMax=" << cFluxXMax
                       << " scalar_flux_net=" << cFluxNet
                       << " scalar_budget_residual=" << cBudgetResidual
                       << Foam::nl;
        }

        if (activityMode == "htHto")
        {
            // A minimal two-species activity transport prototype. HT is
            // advanced first; the HTO source uses HT at the NEW time level,
            // matching implicit HT conversion loss in the discrete budget.
            Foam::volScalarField& HT = htField();
            Foam::volScalarField& HTO = htoField();
            const Foam::scalar oldHT =
                Foam::fvc::domainIntegrate(HT).value();
            const Foam::scalar oldHTO =
                Foam::fvc::domainIntegrate(HTO).value();
            // Integrate a finite-duration emission over the actual Euler
            // time-step overlap, so straddling sourceStart/sourceEnd does not
            // alter the intended total released activity. Uniform physical
            // volume is explicitly a prototype assumption, not stack CFD.
            const Foam::scalar dtSource = runTime.deltaTValue();
            const Foam::scalar tNew = runTime.value();
            const Foam::scalar overlapSeconds = hasActivitySource
              ? Foam::max
                (
                    Foam::scalar(0),
                    Foam::min(tNew, sourceEnd.value())
                  - Foam::max(tNew - dtSource, sourceStart.value())
                )
              : Foam::scalar(0);
            const Foam::scalar activeFraction = overlapSeconds/dtSource;
            const Foam::scalar injectedHT =
                sourceHTRate.value()*overlapSeconds;
            const Foam::scalar injectedHTO =
                sourceHTORate.value()*overlapSeconds;
            const Foam::dimensionSet sourceDensityDims
            (
                0, -3, -2, 0, 0, 0, 0
            );
            const Foam::dimensionedScalar htSourceDensity
            (
                "htSourceDensity", sourceDensityDims,
                sourceHTRate.value()*activeFraction/sourceDomainVolume
            );
            const Foam::dimensionedScalar htoSourceDensity
            (
                "htoSourceDensity", sourceDensityDims,
                sourceHTORate.value()*activeFraction/sourceDomainVolume
            );
            if (hasActivitySource)
            {
                Foam::Info << "CTEAS_M1C emitted_HT_Bq=" << injectedHT
                           << " emitted_HTO_Bq=" << injectedHTO
                           << " source_overlap_s=" << overlapSeconds
                           << Foam::nl;
            }
            // Correct the coupled HT/HTO boundary values (MPI halo exchange)
            // before assembling a corrected non-orthogonal Laplacian. In
            // particular, processorCyclic must use the current cell states
            // when evaluating the explicit snGrad correction. Omitting this
            // yielded partition-dependent species fields even with fully
            // converged linear solves; the default diffusion scheme stays
            // Gauss linear corrected and remains independently configurable.
            HT.correctBoundaryConditions();
            HTO.correctBoundaryConditions();
            // Convexity proof requires the low solution be nonnegative,
            // and both FVM equations share old time, source and reaction.
            const auto blendWithLow = [&](Foam::volScalarField& speciesField,
                const Foam::scalarField& low, const char* species)->Foam::scalar
            {
                const Foam::scalarField& high = speciesField.primitiveField();
                Foam::scalar alpha = 1.0;
                bool restricted = false;
                forAll(low, cellI)
                {
                    const Foam::scalar lo = low[cellI];
                    const Foam::scalar hi = high[cellI];
                    if (!std::isfinite(lo) || lo < 0 || !std::isfinite(hi))
                    {
                        FatalErrorInFunction
                            << "M1C18_GLOBAL_BLEND_INVALID species=" << species
                            << " cell=" << cellI << " low=" << lo
                            << " high=" << hi << exit(FatalError);
                    }
                    if (hi < 0)
                    {
                        restricted = true;
                        alpha = Foam::min(alpha, lo/(lo - hi));
                    }
                }
                // Leave a small margin against rounding down to negative.
                if (restricted) alpha *= (1.0 - 1e-12);
                if (!std::isfinite(alpha) || alpha < 0 || alpha > 1)
                {
                    FatalErrorInFunction
                        << "M1C18_GLOBAL_BLEND_INVALID alpha=" << alpha
                        << exit(FatalError);
                }
                Foam::scalarField& target = speciesField.primitiveFieldRef();
                forAll(target, cellI)
                {
                    target[cellI] = (1.0-alpha)*low[cellI]
                                  + alpha*target[cellI];
                    if (!std::isfinite(target[cellI]) || target[cellI] < -1e-12)
                    {
                        FatalErrorInFunction
                            << "M1C18_GLOBAL_BLEND_INVALID blended species=" << species
                            << " cell=" << cellI << " value=" << target[cellI]
                            << exit(FatalError);
                    }
                }
                speciesField.correctBoundaryConditions();
                Foam::Info << "CTEAS_M1C18 species=" << species
                           << " global_alpha=" << alpha
                           << " limited=" << Foam::label(restricted)
                           << Foam::nl;
                return alpha;
            };
            Foam::autoPtr<Foam::surfaceScalarField> htBlendFlux;
            Foam::autoPtr<Foam::surfaceScalarField> htoBlendFlux;
            Foam::fvScalarMatrix htEqn
            (
                Foam::fvm::ddt(HT)
              + Foam::fvm::div(phi, HT)
              - Foam::fvm::laplacian(DHT, HT)
              + Foam::fvm::Sp(kOx + lambdaT, HT)
             == htSourceDensity*sourceProfile
            );
            htEqn.solve();
            if (activityConvectionBlend == "globalPositive")
            {
                const Foam::scalarField lowHT(HT.primitiveField());
                const Foam::surfaceScalarField lowFlux(htEqn.flux());
                HT.correctBoundaryConditions();
                Foam::fvScalarMatrix htHighEqn
                (
                    Foam::fvm::ddt(HT)
                  + Foam::fvm::div(phi, HT, "div(phi,HTHigh)")
                  - Foam::fvm::laplacian(DHT, HT)
                  + Foam::fvm::Sp(kOx + lambdaT, HT)
                 == htSourceDensity*sourceProfile
                );
                htHighEqn.solve();
                const Foam::surfaceScalarField highFlux(htHighEqn.flux());
                const Foam::scalar alpha = blendWithLow(HT,lowHT,"HT");
                htBlendFlux.reset
                (
                    new Foam::surfaceScalarField
                    (
                        (1.0-alpha)*lowFlux + alpha*highFlux
                    )
                );
            }
            Foam::fvScalarMatrix htoEqn
            (
                Foam::fvm::ddt(HTO)
              + Foam::fvm::div(phi, HTO)
              - Foam::fvm::laplacian(DHTO, HTO)
              + Foam::fvm::Sp(lambdaT, HTO)
             == kOx*HT + htoSourceDensity*sourceProfile
            );
            htoEqn.solve();
            if (activityConvectionBlend == "globalPositive")
            {
                const Foam::scalarField lowHTO(HTO.primitiveField());
                const Foam::surfaceScalarField lowFlux(htoEqn.flux());
                HTO.correctBoundaryConditions();
                Foam::fvScalarMatrix htoHighEqn
                (
                    Foam::fvm::ddt(HTO)
                  + Foam::fvm::div(phi, HTO, "div(phi,HTOHigh)")
                  - Foam::fvm::laplacian(DHTO, HTO)
                  + Foam::fvm::Sp(lambdaT, HTO)
                 == kOx*HT + htoSourceDensity*sourceProfile
                );
                htoHighEqn.solve();
                const Foam::surfaceScalarField highFlux(htoHighEqn.flux());
                const Foam::scalar alpha = blendWithLow(HTO,lowHTO,"HTO");
                htoBlendFlux.reset
                (
                    new Foam::surfaceScalarField
                    (
                        (1.0-alpha)*lowFlux + alpha*highFlux
                    )
                );
            }
            // Fail closed on a genuine unsafe time step BEFORE any output is
            // accepted. Do not alter HT/HTO fields, Bq budgets or fluxes.
            requireNonnegativeActivity("postSolve");
            // M1C-16 opt-in per-step global native primitive-cell
            // extrema. No added flux or transport correction. These
            // values do not imply positivity of interpolated faces.
            if
            (
                transportProperties.lookupOrDefault<Foam::Switch>
                (
                    "auditActivityExtrema", false
                )
            )
            {
                const Foam::scalar htMin =
                    Foam::gMin(HT.primitiveField());
                const Foam::scalar htMax =
                    Foam::gMax(HT.primitiveField());
                const Foam::scalar htoMin =
                    Foam::gMin(HTO.primitiveField());
                const Foam::scalar htoMax =
                    Foam::gMax(HTO.primitiveField());
                // MessageStream does not expose Ostream::precision
                // in official v2606; serialize at full precision without
                // altering formatting of any previous diagnostic log.
                std::ostringstream auditLine;
                auditLine << std::setprecision(16)
                          << "CTEAS_M1C extrema_HT_min_Bq_m3="
                          << htMin
                          << " extrema_HT_max_Bq_m3=" << htMax
                          << " extrema_HTO_min_Bq_m3=" << htoMin
                          << " extrema_HTO_max_Bq_m3=" << htoMax;
                Foam::Info << auditLine.str().c_str() << Foam::nl;
            }
            const Foam::scalar newHT =
                Foam::fvc::domainIntegrate(HT).value();
            const Foam::scalar newHTO =
                Foam::fvc::domainIntegrate(HTO).value();
            Foam::Info << "CTEAS_M1C activity_HT_Bq=" << newHT
                       << " activity_HTO_Bq=" << newHTO << Foam::nl;
            // Optional equation-level activity budget. Use the solved FVM
            // matrix face fluxes, including the selected advective boundary
            // coefficients and diffusion. Requires fluxRequired(HT,HTO).
            // No assertion that a cyclic or open boundary has zero net flux.
            if
            (
                transportProperties.lookupOrDefault<Foam::Switch>
                (
                    "auditActivityFlux", false
                )
            )
            {
                const Foam::label leftPatch =
                    mesh.boundaryMesh().findPatchID("xMin");
                const Foam::label rightPatch =
                    mesh.boundaryMesh().findPatchID("xMax");
                if (leftPatch < 0 || rightPatch < 0)
                {
                    FatalErrorInFunction
                        << "auditActivityFlux requires xMin/xMax patches"
                        << exit(FatalError);
                }
                auto boundaryFlux = [&]
                (
                    const Foam::surfaceScalarField& f,
                    Foam::scalar& left,
                    Foam::scalar& right
                ) -> Foam::scalar
                {
                    Foam::scalar net = 0.0;
                    left = 0.0;
                    right = 0.0;
                    forAll(f.boundaryField(), patchI)
                    {
                        const Foam::scalar value =
                            Foam::sum(f.boundaryField()[patchI]);
                        net += value;
                        // A decomposition crossing a cyclic pair converts
                        // physical xMin/xMax faces into processorCyclic
                        // patches. Their original identity is recorded by
                        // the official polyPatch referPatchName() API.
                        // Count them in the physical signed boundary flux,
                        // not as an ordinary processor internal interface.
                        Foam::word physicalName =
                            mesh.boundaryMesh()[patchI].name();
                        if
                        (
                            Foam::isA<Foam::processorCyclicPolyPatch>
                            (
                                mesh.boundaryMesh()[patchI]
                            )
                        )
                        {
                            physicalName = Foam::refCast
                            <const Foam::processorCyclicPolyPatch>
                            (
                                mesh.boundaryMesh()[patchI]
                            ).referPatchName();
                        }
                        if (physicalName == "xMin") left += value;
                        if (physicalName == "xMax") right += value;
                    }
                    Foam::reduce(net, Foam::sumOp<Foam::scalar>());
                    Foam::reduce(left, Foam::sumOp<Foam::scalar>());
                    Foam::reduce(right, Foam::sumOp<Foam::scalar>());
                    return net;
                };
                Foam::surfaceScalarField htFlux(htEqn.flux());
                Foam::surfaceScalarField htoFlux(htoEqn.flux());
                if (activityConvectionBlend == "globalPositive")
                {
                    htFlux = htBlendFlux();
                    htoFlux = htoBlendFlux();
                }
                Foam::scalar htLeft = 0.0, htRight = 0.0;
                Foam::scalar htoLeft = 0.0, htoRight = 0.0;
                const Foam::scalar htNet =
                    boundaryFlux(htFlux, htLeft, htRight);
                const Foam::scalar htoNet =
                    boundaryFlux(htoFlux, htoLeft, htoRight);
                const Foam::scalar dt = runTime.deltaTValue();
                const Foam::scalar htResidual =
                    newHT-oldHT + dt*htNet
                  + dt*(kOx.value()+lambdaT.value())*newHT
                  - injectedHT;
                const Foam::scalar htoResidual =
                    newHTO-oldHTO + dt*htoNet
                  + dt*lambdaT.value()*newHTO-dt*kOx.value()*newHT
                  - injectedHTO;
                const Foam::scalar totalResidual =
                    newHT+newHTO-oldHT-oldHTO + dt*(htNet+htoNet)
                  + dt*lambdaT.value()*(newHT+newHTO)
                  - injectedHT-injectedHTO;
                Foam::Info << "CTEAS_M1C flux_HT_xMin_Bq_s=" << htLeft
                           << " flux_HT_xMax_Bq_s=" << htRight
                           << " flux_HT_net_Bq_s=" << htNet
                           << " flux_HT_budget_residual_Bq=" << htResidual
                           << Foam::nl;
                Foam::Info << "CTEAS_M1C flux_HTO_xMin_Bq_s=" << htoLeft
                           << " flux_HTO_xMax_Bq_s=" << htoRight
                           << " flux_HTO_net_Bq_s=" << htoNet
                           << " flux_HTO_budget_residual_Bq=" << htoResidual
                           << " flux_total_budget_residual_Bq="
                           << totalResidual << Foam::nl;
            }
            if
            (
                transportProperties.lookupOrDefault<Foam::Switch>
                (
                    "auditClosedActivity", false
                )
            )
            {
                // Valid only for a truly closed, zero-flux box; do not
                // interpret these as open-boundary activity budgets.
                const Foam::scalar dt = runTime.deltaTValue();
                const Foam::scalar rHT =
                    newHT - oldHT + dt*(kOx.value()+lambdaT.value())*newHT
                  - injectedHT;
                const Foam::scalar rHTO =
                    newHTO - oldHTO
                  + dt*lambdaT.value()*newHTO - dt*kOx.value()*newHT
                  - injectedHTO;
                const Foam::scalar rTotal =
                    newHT + newHTO - oldHT - oldHTO
                  + dt*lambdaT.value()*(newHT + newHTO)
                  - injectedHT-injectedHTO;
                Foam::Info << "CTEAS_M1C closed_HT_residual_Bq=" << rHT
                           << " closed_HTO_residual_Bq=" << rHTO
                           << " closed_total_residual_Bq=" << rTotal
                           << Foam::nl;
            }
        }

        Foam::Info << "CTEAS_M1A integral="
                   << Foam::fvc::domainIntegrate(c).value() << Foam::nl;
        runTime.write();
    }

    Foam::Info << "Completed mode=" << momentumMode
               << "; URANS/soil NOT IMPLEMENTED; full HT/HTO physics NOT VALIDATED"
               << Foam::nl;
    return 0;
}
