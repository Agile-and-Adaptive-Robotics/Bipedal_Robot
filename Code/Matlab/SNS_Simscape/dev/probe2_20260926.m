function probe2_20260926()
sns = fileparts(fileparts(mfilename('fullpath')));
mdl = 'mdl_humanoid_lower_ah001_imported';
sub = [mdl '/x09_BA_001_1'];
load_system(fullfile(sns, [mdl '.slx']));
pe = [sub '/x02_01_PE_001_1_RIGID'];
ph = get_param(pe, 'PortHandles');
ports = [ph.RConn ph.LConn];
for pp = 3:8
    porth = ports(pp);
    ln = get_param(porth, 'Line');
    fprintf('port%d line=%.0f', pp, ln);
    try, fprintf(' SrcPH=%.0f', get_param(ln, 'SrcPortHandle')); catch, fprintf(' SrcPH=ERR'); end
    try, fprintf(' DstPH=%.0f', get_param(ln, 'DstPortHandle')); catch, fprintf(' DstPH=ERR'); end
    try, fprintf(' BRH=%s', mat2str(get_param(ln, 'BranchHandles'))); catch, end
    fprintf('\n');
    for fld = ["SrcPortHandle", "DstPortHandle"]
        try
            h = get_param(ln, char(fld));
            if isscalar(h) && h > 0 && h ~= porth
                fprintf('   far via %s: parent=%s\n', char(fld), get_param(get_param(h, 'Parent'), 'Name'));
            end
        catch
        end
    end
end
close_system(mdl, 0);
end
