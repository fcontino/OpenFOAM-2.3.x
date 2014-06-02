tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    reaction->correct();
    dQ = reaction->dQ();
    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        bool tmpActive = composition.isActive(i);
        //if not active on this proc, try to see if it is on another
        //otherwise an error will occur because one of the procs will
        //try to solve the transport equation and the other not
        //if in serial the reduce function return immediately
        reduce(composition.isActive(i), orOp<bool>());

        //composition.isActive(i) returns a reference, here it should be false
        //if all processors have it false or true if at least one is true
        //in the latter case, we should activate the species on the processors
        //where it was inactive
        if (tmpActive != composition.isActive(i))
        {
            composition.setActive(i);
        }

        if (Y[i].name() != inertSpecie && composition.isActive(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->muEff(), Yi)
             ==
                reaction->R(Yi)
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
        else if (Y[i].name() == inertSpecie)
        {
            inertIndex = i;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
