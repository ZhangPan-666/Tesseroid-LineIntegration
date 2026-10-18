function results = test1_run_case(caseName,heightAboveTopMetres)
%TEST1_RUN_CASE Checkpointed 20-component publication calculation.

here = fileparts(mfilename("fullpath"));
kernelDirectory = here;
project = fileparts(fileparts(here));
functionsDirectory = fullfile(project,"Functions");
addpath(kernelDirectory,functionsDirectory);

geometry = test1_geometry(heightAboveTopMetres);
settings = test1_settings();
meta = t1s_components();
caseDirectory = fullfile(here,"results",caseName);
if ~isfolder(caseDirectory)
    mkdir(caseDirectory);
end
checkpointFile = fullfile(caseDirectory,"checkpoint.mat");
progressFile = fullfile(caseDirectory,"progress.log");

if isfile(checkpointFile)
    loaded = load(checkpointFile,"results");
    results = loaded.results;
    local_report(progressFile,"RESUME from checkpoint");
else
    results.schemaVersion = 1;
    results.caseName = caseName;
    results.created = datetime("now");
    results.geometry = geometry;
    results.settings = settings;
    results.components = meta;
    results.environment = local_environment(functionsDirectory, ...
        kernelDirectory);
    results.cpp = struct();
    for methodName = settings.methodOrder
        results.vpa.(methodName) = cell(20,1);
    end
    local_report(progressFile,"NEW run");
    save(checkpointFile,"results","-v7.3");
end

local_report(progressFile,sprintf( ...
    "CASE %s: lat=%.12g deg, lon=%.12g deg, height=%.12g m", ...
    caseName,geometry.point(1),geometry.point(2),heightAboveTopMetres));
local_report(progressFile,sprintf( ...
    "REQUESTED tolerances: AbsTol=%.3e RelTol=%.3e", ...
    settings.requestedAbsTol,settings.requestedRelTol));

if ~isfield(results.cpp,"completed") || ~results.cpp.completed
    local_report(progressFile,"C++ START: parallel entry, all 20 components");
    results.cpp = local_cpp20(geometry,settings,meta);
    results.cpp.completed = true;
    save(checkpointFile,"results","-v7.3");
    local_report(progressFile,sprintf( ...
        "C++ DONE: %.6f s; warning=%s",results.cpp.elapsedSeconds, ...
        results.cpp.warningMessage));
    local_write_tables(results,caseDirectory);
end

methods = {@t1s_vpa1d_component,@t1s_vpa2d_component, ...
    @t1s_vpa3d_component};
totalLevels = 20*numel(methods)*numel(settings.workingDigits);
runStart = tic;
for methodIndex = 1:numel(settings.methodOrder)
    methodName = settings.methodOrder(methodIndex);
    for componentIndex = 1:20
        run = results.vpa.(methodName){componentIndex};
        if ~isstruct(run)
            run = local_empty_run(settings);
        end
        for levelIndex = 1:numel(settings.workingDigits)
            if run.levels(levelIndex).completed
                continue
            end
            completedBefore = local_completed_levels(results);
            elapsedKnown = local_known_vpa_seconds(results);
            eta = local_eta(elapsedKnown,completedBefore,totalLevels);
            wd = settings.workingDigits(levelIndex);
            local_report(progressFile,sprintf( ...
                "START %s %s digits=%d level=%d/%d progress=%d/%d ETA=%s", ...
                methodName,meta.names(componentIndex),wd,levelIndex, ...
                numel(settings.workingDigits),completedBefore,totalLevels,eta));

            started = tic;
            try
                value = methods{methodIndex}(geometry.prism,geometry.point, ...
                    componentIndex,wd,settings.quadratureDigits);
                message = "";
            catch exception
                value = sym(NaN);
                message = string(getReport(exception,"extended", ...
                    "hyperlinks","off"));
            end
            levelSeconds = toc(started);
            run.levels(levelIndex).value = value;
            run.levels(levelIndex).valueText = string(char(vpa(value,wd)));
            run.levels(levelIndex).elapsedSeconds = levelSeconds;
            run.levels(levelIndex).message = message;
            run.levels(levelIndex).finite = ...
                isfinite(double(vpa(value,16)));
            run.levels(levelIndex).completed = true;

            if levelIndex > 1 && run.levels(levelIndex).finite && ...
                    run.levels(levelIndex-1).finite
                prior = run.levels(levelIndex-1).value;
                difference = abs(value-prior);
                relTol = vpa(settings.compareRelTol,wd);
                absTol = vpa(settings.compareAbsTol,wd);
                criterion = max(absTol,relTol*max(abs(value),abs(prior)));
                run.levels(levelIndex).difference = difference;
                run.levels(levelIndex).criterion = criterion;
                run.levels(levelIndex).stable = ...
                    isAlways(difference<=criterion);
            end
            results.vpa.(methodName){componentIndex} = run;
            save(checkpointFile,"results","-v7.3");
            local_write_history(results,caseDirectory);
            local_report(progressFile,sprintf( ...
                "DONE  %s %s digits=%d seconds=%.3f stable=%d message=%s", ...
                methodName,meta.names(componentIndex),wd,levelSeconds, ...
                run.levels(levelIndex).stable,message));
        end

        stable = [run.levels.stable];
        run.converged = all(stable(end-settings.requiredStableSteps+1:end));
        run.value = run.levels(end).value;
        run.valueText = run.levels(end).valueText;
        run.usedDigits = settings.workingDigits(end);
        run.elapsedSeconds = sum([run.levels.elapsedSeconds]);
        run.completed = true;
        results.vpa.(methodName){componentIndex} = run;
        save(checkpointFile,"results","-v7.3");
        local_write_tables(results,caseDirectory);
        local_report(progressFile,sprintf( ...
            "COMPONENT DONE %s %s converged=%d cumulative=%.3f s", ...
            methodName,meta.names(componentIndex),run.converged, ...
            run.elapsedSeconds));
    end
end

results.completed = datetime("now");
results.resumeSessionSeconds = toc(runStart);
results.summary = local_summary(results);
save(checkpointFile,"results","-v7.3");
local_write_tables(results,caseDirectory);
local_report(progressFile,"RUN COMPLETE");
disp(results.summary);
end

function run = local_empty_run(settings)
blank = struct("completed",false,"finite",false,"stable",false, ...
    "value",sym(NaN),"valueText","","difference",sym(NaN), ...
    "criterion",sym(NaN),"elapsedSeconds",NaN,"message","");
run.levels = repmat(blank,numel(settings.workingDigits),1);
run.completed = false;
run.converged = false;
run.value = sym(NaN);
run.valueText = "";
run.usedDigits = NaN;
run.elapsedSeconds = NaN;
end

function out = local_cpp20(geometry,settings,meta)
p = geometry.prism;
model = [p.lat1,p.lat2,p.lon1,p.lon2,p.r1,p.r2,p.density];
point = geometry.point;
lastwarn("");
started = tic;
[V,Vx,Vy,Vz,Vxx,Vxy,Vyy,Vxz,Vyz,Vzz, ...
    Vxxx,Vxxy,Vxxz,Vxyz,Vxyy,Vyyy,Vyyz,Vxzz,Vyzz,Vzzz] = ...
    Tesseroid_GravityEstimateP(model,point, ...
    settings.requestedAbsTol,settings.requestedRelTol);
out.elapsedSeconds = toc(started);
[warningMessage,warningIdentifier] = lastwarn;
out.warningMessage = string(warningMessage);
out.warningIdentifier = string(warningIdentifier);
out.valuesPhysical = [V,Vx,Vy,Vz,Vxx,Vxy,Vxz,Vyy,Vyz,Vzz, ...
    Vxxx,Vxxy,Vxxz,Vxyy,Vxyz,Vxzz,Vyyy,Vyyz,Vyzz,Vzzz];
out.normalizationScale = max(point(3),p.r2);
out.valuesNormalized = out.valuesPhysical./ ...
    out.normalizationScale.^(2-meta.orders.');
out.requestedAbsTol = settings.requestedAbsTol;
out.requestedRelTol = settings.requestedRelTol;
out.effectiveRelTolFloor = settings.cppEffectiveRelTolFloor;
out.entryPoint = "Tesseroid_GravityEstimateP";
end

function environment = local_environment(functionsDirectory,kernelDirectory)
environment.matlabVersion = string(version);
environment.release = string(version("-release"));
environment.computer = string(computer);
environment.numberOfProcessors = string(getenv("NUMBER_OF_PROCESSORS"));
environment.started = datetime("now");
environment.mexPath = string(which("tli_lineintegration_mex"));
environment.mexFile = local_file_info(fullfile(functionsDirectory, ...
    "tli_lineintegration_mex.mexw64"));
environment.estimateHeader = local_file_info(fullfile(fileparts( ...
    functionsDirectory),"C++Codes","_TLI_estimate.h"));
environment.kernel20File = local_file_info(fullfile(kernelDirectory, ...
    "t1s_kernel20.m"));
environment.matlabExecutable = string(matlabroot);
end

function info = local_file_info(path)
d = dir(path);
info.path = string(path);
if isempty(d)
    info.exists = false;
    info.bytes = NaN;
    info.modified = NaT;
else
    info.exists = true;
    info.bytes = d.bytes;
    info.modified = datetime(d.datenum,"ConvertFrom","datenum");
end
end

function n = local_completed_levels(results)
n = 0;
for methodName = results.settings.methodOrder
    for k = 1:20
        run = results.vpa.(methodName){k};
        if isstruct(run) && isfield(run,"levels")
            n = n+sum([run.levels.completed]);
        end
    end
end
end

function seconds = local_known_vpa_seconds(results)
seconds = 0;
for methodName = results.settings.methodOrder
    for k = 1:20
        run = results.vpa.(methodName){k};
        if isstruct(run) && isfield(run,"levels")
            values = [run.levels.elapsedSeconds];
            seconds = seconds+sum(values(isfinite(values)));
        end
    end
end
end

function text = local_eta(elapsed,completed,total)
if completed==0 || elapsed<=0
    text = "unknown";
else
    remaining = elapsed/completed*(total-completed);
    text = string(duration(0,0,remaining,"Format","hh:mm:ss"));
end
end

function local_report(progressFile,message)
line = sprintf("[%s] %s", ...
    char(datetime("now","Format","yyyy-MM-dd HH:mm:ss")),message);
fprintf("%s\n",line);
fid = fopen(progressFile,"a");
if fid>=0
    fprintf(fid,"%s\n",line);
    fclose(fid);
end
end

function local_write_history(results,caseDirectory)
rows = {};
for methodName = results.settings.methodOrder
    for componentIndex = 1:20
        run = results.vpa.(methodName){componentIndex};
        if ~isstruct(run) || ~isfield(run,"levels")
            continue
        end
        for levelIndex = 1:numel(run.levels)
            level = run.levels(levelIndex);
            if ~level.completed
                continue
            end
            rows(end+1,:) = {methodName, ...
                results.components.names(componentIndex),componentIndex, ...
                results.settings.workingDigits(levelIndex), ...
                level.valueText,level.elapsedSeconds,level.finite, ...
                level.stable,string(char(vpa(level.difference,8))), ...
                string(char(vpa(level.criterion,8))),level.message}; %#ok<AGROW>
        end
    end
end
if isempty(rows)
    return
end
history = cell2table(rows,"VariableNames", ...
    ["method","component","componentIndex","workingDigits","valueText", ...
    "elapsedSeconds","finite","stable","difference","criterion","message"]);
writetable(history,fullfile(caseDirectory,"precision_history.csv"));
end

function local_write_tables(results,caseDirectory)
summary = local_summary(results);
writetable(summary,fullfile(caseDirectory,"summary_current.csv"));
metadata = results.environment;
metadata.caseName = results.caseName;
metadata.geometry = results.geometry;
metadata.settings = results.settings;
fid = fopen(fullfile(caseDirectory,"metadata.json"),"w");
if fid>=0
    fprintf(fid,"%s",jsonencode(metadata,"PrettyPrint",true));
    fclose(fid);
end
end

function summary = local_summary(results)
rows = cell(20,25);
for k = 1:20
    one = local_final(results.vpa.gaugeLine1d{k});
    two = local_final(results.vpa.radialSurface2d{k});
    three = local_final(results.vpa.volume3d{k});
    cpp = results.cpp.valuesNormalized(k);
    [a32,r32] = local_error(three.value,two.value);
    [a21,r21] = local_error(two.value,one.value);
    [a31,r31] = local_error(three.value,one.value);
    [a3c,r3c] = local_error(three.value,cpp);
    [~,r1c] = local_error(one.value,cpp);
    rows(k,:) = {results.components.names(k),results.components.orders(k), ...
        one.completed,one.converged,one.valueText,one.seconds, ...
        two.completed,two.converged,two.valueText,two.seconds, ...
        three.completed,three.converged,three.valueText,three.seconds, ...
        results.cpp.valuesPhysical(k),cpp,a32,r32,a21,r21,a31,r31, ...
        a3c,r3c,r1c};
end
summary = cell2table(rows,"VariableNames", ...
    ["component","order","done1d","conv1d","value1d","seconds1d", ...
    "done2d","conv2d","value2d","seconds2d", ...
    "done3d","conv3d","value3d","seconds3d", ...
    "cppPhysical","cppNormalized","abs3d2d","rel3d2d", ...
    "abs2d1d","rel2d1d","abs3d1d","rel3d1d", ...
    "abs3dCpp","rel3dCpp","rel1dCpp"]);
end

function out = local_final(run)
out.completed = false;
out.converged = false;
out.value = sym(NaN);
out.valueText = "";
out.seconds = NaN;
if isstruct(run) && isfield(run,"completed") && run.completed
    out.completed = true;
    out.converged = run.converged;
    out.value = run.value;
    out.valueText = run.valueText;
    out.seconds = run.elapsedSeconds;
end
end

function [absolute,relative] = local_error(a,b)
if ~isfinite(double(vpa(a,16))) || ~isfinite(double(vpa(b,16)))
    absolute = NaN;
    relative = NaN;
    return
end
aa = vpa(a,32);
bb = vpa(b,32);
absolute = double(vpa(abs(aa-bb),16));
relative = double(vpa(abs(aa-bb)/max([abs(aa),abs(bb),sym(1e-100)]),16));
end
