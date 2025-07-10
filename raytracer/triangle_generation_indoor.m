function [s] = triangle_generation_indoor(s,idx)
% This function generates triangles on each surfaces in the scenario for 
% reflection plane selection. Each surface is assumed to be square and 
% consists of two triangles. For oindoor scenarios, four walls, a ceiling 
% and a floor should be considered along with the surfaces of scatterers.
% The triangles are described by their vertices.

% Input:
% s: the floor map which has been converted to struct format
% idx: index

% Output:
% s: the floor map

% Author: Yunqi Feng <Yunqi.Feng@UGent.be>
% September 2023
s.amf.object{1,6+idx}.mesh.volume.triangle{1,1}.v1.Text='0';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,1}.v2.Text='1';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,1}.v3.Text='3';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,2}.v1.Text='1';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,2}.v2.Text='3';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,2}.v3.Text='2';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,3}.v1.Text='1';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,3}.v2.Text='2';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,3}.v3.Text='5';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,4}.v1.Text='2';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,4}.v2.Text='5';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,4}.v3.Text='6';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,5}.v1.Text='4';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,5}.v2.Text='5';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,5}.v3.Text='6';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,6}.v1.Text='4';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,6}.v2.Text='6';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,6}.v3.Text='7';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,7}.v1.Text='0';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,7}.v2.Text='3';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,7}.v3.Text='4';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,8}.v1.Text='3';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,8}.v2.Text='4';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,8}.v3.Text='7';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,9}.v1.Text='2';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,9}.v2.Text='3';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,9}.v3.Text='6';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,10}.v1.Text='3';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,10}.v2.Text='6';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,10}.v3.Text='7';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,11}.v1.Text='0';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,11}.v2.Text='1';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,11}.v3.Text='5';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,12}.v1.Text='0';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,12}.v2.Text='4';
s.amf.object{1,6+idx}.mesh.volume.triangle{1,12}.v3.Text='5';
end