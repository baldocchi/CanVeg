function [quantum, nir, ir, Qin, rnet, Sun, Shade,Veg]=fZeroArrays(prm)

% zero arrays for speed     
%  Dec 3, 2020

      quantum.sh_abs=zeros(prm.nn,prm.jtot);
      quantum.sun_abs=zeros(prm.nn,prm.jtot);
      
      quantum.sun=zeros(prm.nn,prm.jktot);
      quantum.dn_flux=zeros(prm.nn,prm.jktot);
      quantum.up_flux=zeros(prm.nn,prm.jktot);
      quantum.sun_normal=zeros(prm.nn,prm.jktot);
      quantum.sh_flux=zeros(prm.nn,prm.jktot);
      quantum.incoming=zeros(prm.nn,1);
      quantum.beam_flux=zeros(prm.nn,prm.jktot);
      quantum.total=zeros(prm.nn,prm.jktot);
      quantum.inbeam=zeros(prm.nn,1);
      quantum.indiffuse=zeros(prm.nn,1);
      quantum.prob_beam=zeros(prm.nn,1);
      

      
      quantum.prob_beam=zeros(prm.nn,prm.jktot);
      quantum.prob_shade=ones(prm.nn,prm.jktot);      % initialize prob of shade to one if night
      quantum.sun_lai=zeros(prm.nn,prm.jktot);
      quantum.shade_lai=zeros(prm.nn,prm.jktot);
      
      
%       par.sh_abs=zeros(prm.nn,prm.jtot);
%       par.sun_abs=zeros(prm.nn,prm.jtot);
%       
%       par.sun=zeros(prm.nn,prm.jktot);
%       par.dn_flux=zeros(prm.nn,prm.jktot);
%       par.up_flux=zeros(prm.nn,prm.jktot);
%       par.sun_normal=zeros(prm.nn,prm.jktot);
%       par.sh_flux=zeros(prm.nn,prm.jktot);
%       par.incoming=zeros(prm.nn,1);
%       par.beam_flux=zeros(prm.nn,prm.jktot);
%       par.total=zeros(prm.nn,prm.jktot);
%       
%       par.sun_top_flux=zeros(prm.nn,prm.jktot);
%       par.sh_top_flux=zeros(prm.nn,prm.jktot);
%       par.sh_bottom_flux=zeros(prm.nn,prm.jktot);
%       
%       par.prob_beam=zeros(prm.nn,prm.jktot);
%       par.prob_shade=ones(prm.nn,prm.jktot);      % initialize prob of shade to one if night
%       par.sun_lai=zeros(prm.nn,prm.jktot);
%       par.shade_lai=zeros(prm.nn,prm.jktot);
      
      nir.sh_abs=zeros(prm.nn,prm.jtot);
      nir.sun_abs=zeros(prm.nn,prm.jtot);
      
      nir.sun=zeros(prm.nn,prm.jktot);
      nir.dn_flux=zeros(prm.nn,prm.jktot);
      nir.up_flux=zeros(prm.nn,prm.jktot);
      nir.sun_normal=zeros(prm.nn,prm.jktot);
      nir.sh_flux=zeros(prm.nn,prm.jktot);
      
      
      nir.incoming=zeros(prm.nn,1);
      nir.beam_flux=zeros(prm.nn,prm.jktot);
      nir.total=zeros(prm.nn,prm.jktot);
      nir.inbeam=zeros(prm.nn,1);
      nir.indiffuse=zeros(prm.nn,1);
      
      nir.prob_beam=zeros(prm.nn,prm.jktot);
      nir.prob_shade=ones(prm.nn,prm.jktot);      % initialize prob of shade to one if night
      nir.sun_lai=zeros(prm.nn,prm.jktot);
      nir.shade_lai=zeros(prm.nn,prm.jktot);
      
      
      ir.ir_dn=ones(prm.nn,prm.jktot);
      ir.ir_up=ones(prm.nn,prm.jktot);
      
      ir.IR_source_sun=zeros(prm.nn,prm.jktot);
      ir.IR_source_shade=zeros(prm.nn,prm.jktot);
      ir.IR_source=zeros(prm.nn,prm.jktot);
      
      ir.shade=zeros(prm.nn,prm.jtot);
      
      ir.shade_top=zeros(prm.nn,prm.jktot);
      ir.shade_bottom=zeros(prm.nn,prm.jktot);
      
      rnet.sun=zeros(prm.nn,prm.jktot);
      rnet.sh=zeros(prm.nn,prm.jktot);
      
      rnet.sun_top=zeros(prm.nn,prm.jktot);
      rnet.sh_top=zeros(prm.nn,prm.jktot);
      rnet.sh_bottom=zeros(prm.nn,prm.jktot);
      
      
      Sun.Ps=zeros(prm.nn,prm.jtot);
      Sun.Resp=zeros(prm.nn,prm.jtot);
      Sun.gs=zeros(prm.nn,prm.jtot);
      Sun.LE=zeros(prm.nn,prm.jtot);
      Sun.H=zeros(prm.nn,prm.jtot);
      Sun.Tsfc=ones(prm.nn,prm.jtot);
      Sun.Tsfc_old=ones(prm.nn,prm.jtot);
      Sun.Tsfc_new=ones(prm.nn,prm.jtot);
      
      Shade.Ps=zeros(prm.nn,prm.jtot);
      Shade.Resp=zeros(prm.nn,prm.jtot);
      Shade.gs=zeros(prm.nn,prm.jtot);
      Shade.LE=zeros(prm.nn,prm.jtot);
      Shade.H=zeros(prm.nn,prm.jtot);
      Shade.Tsfc=ones(prm.nn,prm.jtot);
      Shade.Tsfc_old=ones(prm.nn,prm.jtot);
      Shade.Tsfc_new=ones(prm.nn,prm.jtot);
      
      boundary_layer_res.heat=zeros(prm.nn,prm.jtot);
      boundary_layer_res.vapor=zeros(prm.nn,prm.jtot);
      boundary_layer_res.co2=zeros(prm.nn,prm.jtot);
      
      Qin.sun_abs=zeros(prm.nn,prm.jtot);
      Qin.shade_abs=zeros(prm.nn,prm.jtot);
      
      Veg.Ps=zeros(prm.nn,1);
      Veg.Rd=zeros(prm.nn,1);
      Veg.H=zeros(prm.nn,1);
      Veg.LE=zeros(prm.nn,1);
      Veg.Tsfc=zeros(prm.nn,1);


end