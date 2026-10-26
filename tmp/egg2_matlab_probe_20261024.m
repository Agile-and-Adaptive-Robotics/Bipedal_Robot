% egg2 MATLAB capability probe 2026-10-24: which Simulink/Simscape
% products are installed here. Output only; no models touched.
disp(['HOST: ' getenv('COMPUTERNAME')]);
disp(['VERSION: ' version]);
v = ver;
disp('--- Simulink/Simscape products ---');
for k = 1:numel(v)
    nm = v(k).Name;
    a = contains(nm, 'Simulink', 'IgnoreCase', true);
    b = contains(nm, 'Simscape', 'IgnoreCase', true);
    if a || b
        disp([nm ' | ' v(k).Version]);
    end
end
disp('--- probe done ---');
