function settings = test1_settings()
%TEST1_SETTINGS Publication run configuration.

settings.requestedAbsTol = 1e-16;
settings.requestedRelTol = 1e-16;
settings.workingDigits = [28,40,52];
settings.quadratureDigits = 16;
settings.compareAbsTol = "1e-16";
settings.compareRelTol = "1e-16";
settings.requiredStableSteps = 2;
settings.cppEffectiveRelTolFloor = 500*eps;
settings.methodOrder = ["gaugeLine1d","radialSurface2d","volume3d"];
settings.checkpointAfterEveryPrecisionLevel = true;
end
