fvScalarMatrix pEqn
(
fvm::laplacian(p_)
);
pEqn.setReference(pRefCell, pRefValue);
