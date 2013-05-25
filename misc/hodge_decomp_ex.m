clear
path(path, './toolbox_diffc/toolbox_diffc');
path(path, './toolbox_diffc/toolbox_diffc/toolbox');

rng('default')
n = 100;
options.bound = 'per';
v = randn(n,n,2);
v = perform_vf_normalization( perform_blurring(v,100) );
[v_nocurl,v_nodiv] = compute_hodge_decompositon(v,options);
options.subsampling = 4;
options.display_streamlines = 0;
options.normalize_flow = 1;
options.lgd = { 'Original=Irrot+Incomp' 'Irrotational' 'Incompressible'  };
figure(1)
plot_vf( {v v_nocurl v_nodiv}, [], options );

[v_nocurl2,v_nodiv2] = compute_hodge_decompositon(v_nocurl,options);
options.subsampling = 4;
options.display_streamlines = 0;
options.normalize_flow = 1;
options.lgd = { 'Irrotational Field' 'Irrotational' 'Incompressible'  };
figure(2)
plot_vf( {v_nocurl v_nocurl2 v_nodiv2}, [], options );

[v_nocurl3,v_nodiv3] = compute_hodge_decompositon(v_nodiv,options);
options.subsampling = 4;
options.display_streamlines = 0;
options.normalize_flow = 1;
options.lgd = { 'Incompressible Field' 'Irrotational' 'Incompressible'  };
figure(3)
plot_vf( {v_nodiv v_nocurl3 v_nodiv3}, [], options );
