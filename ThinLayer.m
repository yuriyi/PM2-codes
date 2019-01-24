(* ::Package:: *)

(* ::Input::Initialization:: *)
BeginPackage["ThinLayer`"]


(* ::Input::Initialization:: *)
(* Usage Messaging *)
tlSetLayer::usage="tlSetLayer[layerIndex, {c11, c13, c33, c44, \[Rho]}] creates a 'tlLayer[layerIndex]' object instance with assigned numerical values for cijs.";
tlGetStiffness::usage="tlGetStiffness[layerIndex] outputs a vector of cijs for a given layer.";
tlCijFromThomsen::usage="tlCijFromThomsen[vp0,vs0,\[Epsilon],\[Delta],\[Rho]] outputs {c11, c13, c33, c44, \[Rho]} for a set of Thomsen parameters.";
tlVTISquirtModuli::usage="tlVTISquirtModuli[\[Epsilon],\[Epsilon]f,\[Alpha],\[Tau]0, \[Omega], lenrat][\[Lambda],\[Mu], \[Phi], Kf] returns {c11, c13, c33, c44} for a given frequency.";
tlQLmatrices::usage="tlQmatices[i, p, \[Omega]:0] returns the numbers qa, qb and linear transformations L1, L2 decomposing the wave field of horizontal slowness p, in a layer with index i. The output is in the form {{qa, qb}, L1, L2}";
tlQLmatricesRadicals::usage="[index_Integer,p_,\[Omega]_:1]";
tlLayerReflectionTransmission::usage="tlLayerReflectionTransmission[i, j, p] calculates transmission T and reflection R matrices for an interface between pre-defined layers i and j for horizontal slowness p. It can also be called as tlLayerReflectionTransmission[i, p] in which case an interface between layers i, i+1 is assumed. The output is of the form {T, R}."; 
tlReflectivity::usage="tlReflectivity[{i,j, ...},{di, dj, ...}, slowness, angularFreq] calculates the reflection response matrix RR of a sequence of pre-initiated layers labelled {i,j,...} with thicknesses {di, dj, ...} where the first and last layers are halfspaces and the first and last sequence is ignored.";
tlResponse::usage="tlResponse[{i,j, ...},{di, dj, ...}, slowness, angularFreq,zs,zb] calculates the stress-displacement vector for a stack of layers bounded by two half-spaces. source and receiver are located in the top half-space: source depth: zs, receiver depth: 0, first boundary: zb. the vetor returned is in frequency-slowness domain: {Uz, Srz, Szz, Ur}";
tlResponseVforce::usage="tlResponseVforce[{i,j, ...},{di, dj, ...}, slowness, angularFreq,zs,zb] calculates the stress-displacement vector for a stack of layers bounded by two half-spaces. source and receiver are located in the top half-space: source depth: zs, receiver depth: 0, first boundary: zb. the vetor returned is in frequency-slowness domain: {Uz, Srz, Szz, Ur}. Source type - vertical force";
tlDHTR::usage="tlDHTR[u,\[Nu],r,rmax,\[Omega],BesselVector] returns HT of order \[Nu] of function u[p] at radius r for given frequency \[Omega]. BesselVector is the vector of zeroes of the Bessel funciton, rmax is the maximum radius after which the function is assumed 0. check the definition to see apply apodizing windows";
tlPlotSeismic::usage="[label,field(\[Omega]-x),\[Omega],offsets,ampl,\[Omega]I,model,dlist,{z_source,z_boundary},output_file:0,fRicker:25,maxt,traveltime or a coded event (1P2PSPS2P) ]";
tlPlotDiff::usage="[delta_,offsets_,t_,ampl_,pr_]";
tlPlotTrace::usage="[u_,\[Omega]_,\[Omega]I_,fRicker_:25]";
tlTravelTimes::usage="tlTravelTimes[{i,j, ...},{di, dj, ...}, slowness, zs: source depth, zb: depth of the first boundary] for a given slowness returns pairs of offset-traveltime {xi, t(xi)} of all the primary relfections of reflected and converted waves arranged as { {PP, SS, PS, SP}, {Direct} }. Direct wave traveltimes are for P and S modes.";
tlVisualiseModel::usage="tlVisualiseModel[{i,j, ...},{di, dj, ...}, nSamples] outputs a graphics object visualising different layers with different colours in the time domain.";(*experimental function. Aiming to display layered model with colours for different layers*)
(* Error Reporting tlSetLayer *)
tlSetLayer::nonum="Expecting an array of 4 complex and 1 real numeric quantities in arg `2`";
tlSetLayer::noint="Argument `1` refers to layer label. Positive integer expected";
tlSetLayer::noagrs="tlSetLayer called with wrong number of arguments. tlSetLayer[layerIndex, {c11, c13, c33, c44, \[Rho]}] expected";
(* Error Reporting tlQLmatrices *)
tlQLmatrices::nodef="Layer `1` not defined. Define elastic constants and density for layer `1` using tlSetLayer";
tlQLmatrices::noint="Argument `1` refers to layer label. Integer expected";
(* Error Reporting tlLayerReflectionTransmission*)
tlLayerReflectionTransmission::nodef="Layers `1` and `2` must be defined with tlSetLayer. Define elastic constants and density for layers `1`, `2` using tlSetLayer";
tlLayerReflectionTransmission::noint="Argument `1` refers to layer label. Integer expected";
(* Error Reporting tlReflectivity*)
tlReflectivity::nodef="Layer `1` not defined. Define elastic constants and density for layer `1` using tlSetLayer";
tlReflectivity::nodef="Argument `1` refers to layer label. Integer expected";


(* ::Input::Initialization:: *)
Begin["Private`"]


(* ::Input::Initialization:: *)
tlSetLayer[index_Integer,{c11_?NumericQ, c13_?NumericQ, c33_?NumericQ, c44_?NumericQ, \[Rho]_Real}]/;(Positive[index]):=
 tlLayer[index]=N@{c11,c13,c33,c44,\[Rho]};
tlSetLayer[index_Integer,a:({_,_,_,_,_})]/;(Positive[index]):=(Message[tlSetLayer::nonum,index];$Failed);
tlSetLayer[index_,a:{_,_,_,_,_}]:=(Message[tlSetLayer::noint,index];$Failed);
tlSetLayer[___]:=(Message[tlSetLayer::noargs];$Failed);


(* ::Input::Initialization:: *)
SetAttributes[tlGetStiffness,Listable];
tlGetStiffness[index_Integer]/;ListQ[tlLayer[index]]:=tlLayer[index]


(* ::Input::Initialization:: *)
tlVTISquirtModuli[\[Epsilon]_,\[Epsilon]f_,\[Alpha]_,\[Tau]0_, \[Omega]_, lenrat_][\[Lambda]_,\[Mu]_,\[Rho]_,\[Phi]_, Kf_]:=Module[
{\[Phi]c, \[Phi]f, \[Phi]p, \[Gamma], \[Gamma]2, \[Iota], \[Beta], \[Sigma]c, \[Nu], Kp, Kc, D1, D2, G1, G2, G3, F1, F2,k, L2,L4, \[Tau]f, \[Tau]m,c},
\[Phi]c=(4 \[Pi] \[Epsilon] \[Alpha])/3.;
\[Phi]f=(4 \[Pi] \[Epsilon]f \[Alpha])/3.;
\[Phi]p=\[Phi]-\[Phi]c-\[Phi]f;
Kp=(4.\[Mu])/(3 Kf);
Kc=\[Sigma]c/Kf;
\[Nu]=\[Lambda]/(2(\[Lambda]+\[Mu]));(* Poisson's ratio *)
\[Sigma]c=(\[Pi] \[Mu] \[Alpha])/(2.(1-\[Nu])); (* Crack stiffness *)

\[Phi]p=\[Phi]-\[Phi]c-\[Phi]f;
(*
PRESSURE CALCULATION
*)
\[Gamma]=(3\[Phi]p \[Sigma]c (1+Kp))/(4\[Phi]c \[Mu] (1+Kc));
\[Gamma]2=\[Gamma] (1-\[Nu])/((1+\[Nu])(1+Kp));
\[Iota]=\[Phi]c/(\[Phi]c+\[Alpha] \[Phi]p);
\[Beta]=\[Iota] \[Phi]f/\[Phi]c;
\[Tau]f=lenrat \[Tau]0;
\[Tau]m=\[Tau]0;
D1=(\[Iota] (-I+\[Tau]f \[Omega]-\[Beta] \[Tau]m \[Omega])+3 (1+Kc) \[Gamma]2 ((1+I \[Tau]f \[Omega]) (-I+\[Tau]m \[Omega])+\[Iota] (I-\[Tau]f \[Omega]+\[Beta] \[Tau]m \[Omega])))/(3 (1+Kc) (\[Beta] (-I+(1+(-1+\[Gamma]) \[Iota]) \[Tau]m \[Omega])-(-I+\[Tau]f \[Omega]) (-\[Iota]+\[Gamma] (-1+\[Iota]-I \[Tau]m \[Omega]))));
D2=(\[Beta] (-I+\[Tau]m \[Omega]))/((1+Kc) (\[Beta] (-I+(1+(-1+\[Gamma]) \[Iota]) \[Tau]m \[Omega])-(-I+\[Tau]f \[Omega]) (-\[Iota]+\[Gamma] (-1+\[Iota]-I \[Tau]m \[Omega]))));
G1=(\[Tau]m \[Omega])/((1+Kc) (-I+\[Tau]m \[Omega]));
G2=(I (I-\[Tau]f \[Omega]+\[Beta] \[Tau]m \[Omega]) (\[Iota]+I \[Gamma] \[Iota] \[Tau]m \[Omega]-3 I (1+Kc) \[Gamma]2 (-1+\[Iota]) (-I+\[Tau]m \[Omega])))/(3 (1+Kc) (-I+\[Tau]m \[Omega]) (\[Beta] (-I+(1+(-1+\[Gamma]) \[Iota]) \[Tau]m \[Omega])-(-I+\[Tau]f \[Omega]) (-\[Iota]+\[Gamma] (-1+\[Iota]-I \[Tau]m \[Omega]))));
G3=(\[Beta] (-I+\[Gamma] \[Tau]m \[Omega]))/((1+Kc) (\[Beta] (-I+(1+(-1+\[Gamma]) \[Iota]) \[Tau]m \[Omega])-(-I+\[Tau]f \[Omega]) (-\[Iota]+\[Gamma] (-1+\[Iota]-I \[Tau]m \[Omega]))));
F1=(-3 (1+Kc) \[Gamma]2 (-1+\[Iota]) (-I+\[Tau]m \[Omega])+\[Iota] (-I+\[Gamma] \[Tau]m \[Omega]))/(3 (1+Kc) (\[Beta] (-I+(1+(-1+\[Gamma]) \[Iota]) \[Tau]m \[Omega])-(-I+\[Tau]f \[Omega]) (-\[Iota]+\[Gamma] (-1+\[Iota]-I \[Tau]m \[Omega]))));
F2=(\[Tau]f \[Omega] (\[Gamma]+\[Iota]-\[Gamma] \[Iota]+I \[Gamma] \[Tau]m \[Omega])+\[Beta] (-I+(1+(-1+\[Gamma]) \[Iota]) \[Tau]m \[Omega]))/((1+Kc) (\[Beta] (-I+(1+(-1+\[Gamma]) \[Iota]) \[Tau]m \[Omega])-(-I+\[Tau]f \[Omega]) (-\[Iota]+\[Gamma] (-1+\[Iota]-I \[Tau]m \[Omega]))));
k=\[Lambda] +2/3\[Mu];
L2=k^2+(16 \[Mu]^2)/45;
L4=k^2-(8 \[Mu]^2)/45;
c[11]=\[Lambda]+2\[Mu]
-\[Phi]c(L2/\[Sigma]c+32/15 (1-\[Nu])/(2-\[Nu]) \[Mu]/(\[Pi] \[Alpha])-(L2/\[Sigma]c+k)G1-((3k^2)/\[Sigma]c+3k)G2-((\[Lambda] k)/\[Sigma]c+\[Lambda])G3)
-\[Phi]p(3/(4\[Mu]) (1-\[Nu])/(1+\[Nu]) (3\[Lambda]^2+4\[Lambda] \[Mu]+(36+20\[Nu])/(7-5\[Nu]) \[Mu]^2)-(1+(3k)/(4\[Mu]))(3k D1+\[Lambda] D2))
-\[Phi]f(\[Lambda]^2/\[Sigma]c-3k(\[Lambda]/\[Sigma]c+1)F1-\[Lambda](\[Lambda]/\[Sigma]c+1)F2);
c[33]=\[Lambda]+2\[Mu]
-\[Phi]c(L2/\[Sigma]c+32/15 (1-\[Nu])/(2-\[Nu]) \[Mu]/(\[Pi] \[Alpha])-(L2/\[Sigma]c+k)G1-((3k^2)/\[Sigma]c+3k)G2-((\[Lambda]+2\[Mu] )/\[Sigma]c k+\[Lambda]+2\[Mu])G3)
-\[Phi]p(3/(4\[Mu]) (1-\[Nu])/(1+\[Nu]) (3\[Lambda]^2+4\[Lambda] \[Mu]+(36+20\[Nu])/(7-5\[Nu]) \[Mu]^2)-(1+(3k)/(4\[Mu]))(3k D1+(\[Lambda]+2\[Mu]) D2))
-\[Phi]f((\[Lambda]+2\[Mu])^2/\[Sigma]c-3k((\[Lambda]+2\[Mu])/\[Sigma]c+1)F1-(\[Lambda]+2\[Mu])((\[Lambda]+2\[Mu])/\[Sigma]c+1)F2);
c[44]=\[Mu]
-\[Phi]c(4/15 \[Mu]^2/\[Sigma]c (1-G1)+8/5 (1-\[Nu])/(2-\[Nu]) \[Mu]/(\[Pi] \[Alpha]))
-\[Phi]p(15\[Mu] (1-\[Nu])/(7-5\[Nu]))
-\[Phi]f(4 (1-\[Nu])/(2-\[Nu]) \[Mu]/(\[Pi] \[Alpha]));
c[12]=\[Lambda]
-\[Phi]c(L4/\[Sigma]c-16/15 (1-\[Nu])/(2-\[Nu]) \[Mu]/(\[Pi] \[Alpha])-(L4/\[Sigma]c+k)G1-((3k^2)/\[Sigma]c+3k)G2-((\[Lambda] k)/\[Sigma]c+\[Lambda])G3)
-\[Phi]p(3/(4\[Mu]) (1-\[Nu])/(1+\[Nu]) (3\[Lambda]^2+4\[Lambda] \[Mu]-(4(1+5\[Nu]))/(7-5\[Nu]) \[Mu]^2)-(1+(3k)/(4\[Mu]))(3k D1+\[Lambda] D2))
-\[Phi]f(\[Lambda]^2/\[Sigma]c-3k(\[Lambda]/\[Sigma]c+1)F1-\[Lambda](\[Lambda]/\[Sigma]c+1)F2);
c[13]=\[Lambda]
-\[Phi]c(L4/\[Sigma]c-16/15 (1-\[Nu])/(2-\[Nu]) \[Mu]/(\[Pi] \[Alpha])-(L4/\[Sigma]c+k)G1-((3k^2)/\[Sigma]c+3k)G2-((\[Lambda]+\[Mu] )/\[Sigma]c k+\[Lambda]+\[Mu])G3)
-\[Phi]p(3/(4\[Mu]) (1-\[Nu])/(1+\[Nu]) (3\[Lambda]^2+4\[Lambda] \[Mu]-(4(1+5\[Nu]))/(7-5\[Nu]) \[Mu]^2)-(1+(3k)/(4\[Mu]))(3k D1+(\[Lambda]+\[Mu]) D2))
-\[Phi]f((\[Lambda](\[Lambda]+2\[Mu]))/\[Sigma]c-3k((\[Lambda]+\[Mu])/\[Sigma]c+1)F1-((\[Lambda](\[Lambda]+2\[Mu]))/\[Sigma]c+\[Lambda]+\[Mu])F2);
c[66]=(c[11]-c[12])/2;
{c[11],c[13],c[33],c[44],\[Rho]}//N
]


(* ::Input::Initialization:: *)
tlCijFromThomsen[vp0_,vs0_,\[Epsilon]_,\[Delta]_,\[Rho]_]:=
Module[{c11,c13,c33,c44},
c11=0.I+(1+2\[Epsilon])vp0^2 \[Rho];
c33=0.I+vp0^2 \[Rho];
c13=0.I+(-vs0^2+Sqrt[(vp0-vs0) (vp0+vs0) (-vs0^2+vp0^2 (1+2 \[Delta]))]) \[Rho];
c44=0.I+vs0^2 \[Rho];
{c11,c13,c33,c44,\[Rho]}//N
]


(* ::Input::Initialization:: *)
SetAttributes[tlQLmatrices,Listable];
tlQLmatrices[index_Integer,p_,\[Omega]_:1]/;ListQ[tlLayer[index]]:=
Block[{c, d},
{c[11],c[13],c[33],c[44],d}=tlLayer[index];
Module[
{sgn,imsign,sqrt4,L1,L2,d1,d2,d3,d4,d5,q\[Alpha],q\[Beta],\[Alpha]0,\[Beta]0,\[Eta],\[Delta],\[Sigma]0,S\[Alpha],S\[Beta],R,R1,R2,req1,req2,imq1,imq2,ds,q\[Alpha]s,q\[Beta]s,q\[Alpha]2,q\[Beta]2},
sgn=2(UnitStep[#]-1/2)&;

imsign=Re@#1+I #2 Block[{im=Im@#1},Sign[im]im]&;
sqrt4=Block[{q=#^(1/4)},Re@q+I Im@q]&;

\[Alpha]0=Sqrt[c[33]/d];
\[Beta]0=Sqrt[c[44]/d];
\[Sigma]0=1-c[44]/c[33];
\[Delta]=(c[13]-c[33]+2c[44])/c[33];
\[Eta]=(c[11] c[33]-(c[13]+2c[44])^2)/(2c[33]^2);
R1=2(1-p^2 \[Beta]0^2)(\[Delta]+2p^2 \[Alpha]0^2 \[Eta])^2;
R2=\[Sigma]0+2p^2 \[Beta]0^2 \[Delta]-2p^2 \[Alpha]0^2 (1-2p^2 \[Beta]0^2)\[Eta];
R=R1/(R2+Sqrt[R2^2+2p^2 \[Beta]0^2 R1]);
S\[Alpha]=2\[Delta]+2p^2 \[Alpha]0^2 \[Eta]+R;
S\[Beta]=2(1-p^2 \[Beta]0^2) \[Alpha]0^2/\[Beta]0^2 \[Eta]-R;

q\[Alpha]2=1/\[Alpha]0^2-p^2-p^2 S\[Alpha];
q\[Beta]2=1/\[Beta]0^2-p^2-p^2 S\[Beta];

q\[Alpha]=imsign[Sqrt[q\[Alpha]2],1];
q\[Beta]=imsign[Sqrt[q\[Beta]2],1];

d2=Sqrt[(\[Sigma]0+\[Delta])/(\[Sigma]0+S\[Alpha])];
d3=2\[Beta]0^2 (\[Sigma]0+1/2 (S\[Alpha]+\[Delta]))/(\[Sigma]0+\[Delta]);
d4=Sqrt[(\[Sigma]0-p^2 \[Beta]0^2 (\[Sigma]0+S\[Beta]))/((1-p^2 \[Beta]0^2 (1+S\[Beta]))(\[Sigma]0+\[Delta]))];
d5=(\[Sigma]0-2p^2 \[Beta]0^2 (\[Sigma]0+1/2 (S\[Beta]+\[Delta])))/(\[Sigma]0+\[Delta]);
d1=1/Sqrt[p^2 d3+d5];
ds=Sqrt[d];

{q\[Alpha]s,q\[Beta]s}={sqrt4[q\[Alpha]2],sqrt4[q\[Beta]2]};


L1=d1{
{d2 q\[Alpha]s/ds,1/d4 p/ds/q\[Beta]s},
{d3 d2 p ds q\[Alpha]s,-d5/d4 ds/q\[Beta]s}
};

L2=d1{
{d5/d2 ds/q\[Alpha]s,d3 d4 p ds q\[Beta]s},
{1/d2 p/ds/q\[Alpha]s,-d4 q\[Beta]s/ds}
};
Developer`ToPackedArray/@{{q\[Alpha],q\[Beta]},L1,L2}
]
];
tlQLmatrices[index_Integer,p_]:=(Message[tlQLmatrices::nodef,index];
$Failed);
tlQLmatrices[index_,p_]:=(Message[tlQLmatrices::noint,index];
$Failed);


(* ::Input::Initialization:: *)
(*
<<<<<<<<<<<<<<<<<<< UNDER DEVELOPMENT >>>>>>>>>>>>>>>>>>>>>

SetAttributes[tlQLmatricesRadicals,Listable];
tlQLmatricesRadicals[index_Integer,p_,\[Omega]_:1]/;ListQ[tlLayer[index]]:=
Block[{c, d},
{c[11],c[13],c[33],c[44],d}=tlLayer[index];
Module[
{L1,L2,d1,d2,d3,d4,d5,q\[Alpha],q\[Beta],\[Alpha]0,\[Beta]0,\[Eta],\[Delta],\[Sigma]0,S\[Alpha],S\[Beta],R,R1,R2,req1,req2,imq1,imq2,ds,q\[Alpha]s,q\[Beta]s,q\[Alpha]2,q\[Beta]2,Q},

\[Alpha]0=Sqrt[c[33]/d];
\[Beta]0=Sqrt[c[44]/d];
\[Sigma]0=1-c[44]/c[33];
\[Delta]=(c[13]-c[33]+2c[44])/c[33];
\[Eta]=(c[11] c[33]-(c[13]+2c[44])^2)/(2c[33]^2);
R1=2(1-p^2 \[Beta]0^2)(\[Delta]+2p^2 \[Alpha]0^2 \[Eta])^2;
R2=\[Sigma]0+2p^2 \[Beta]0^2 \[Delta]-2p^2 \[Alpha]0^2 (1-2p^2 \[Beta]0^2)\[Eta];
R=R1/(R2+Sqrt[R2^2+2p^2 \[Beta]0^2 R1]);
S\[Alpha]=2\[Delta]+2p^2 \[Alpha]0^2 \[Eta]+R;
S\[Beta]=2(1-p^2 \[Beta]0^2) \[Alpha]0^2/\[Beta]0^2 \[Eta]-R;

q\[Alpha]2=1/\[Alpha]0^2-p^2-p^2 S\[Alpha];
q\[Beta]2=1/\[Beta]0^2-p^2-p^2 S\[Beta];

Q[q2_,pc2_]:=Module[{a=Re@q2,b=Im@q2},
Which[
p^2<pc2,
{Sqrt[Sqrt[a^2+b^2]+a]+I Sign[\[Omega]]Sqrt[Sqrt[a^2+b^2]-a],
Sqrt[(a^2+b^2)^(1/4)+Sqrt[(Sqrt[a^2+b^2]+a)/2.]]+I Sign[\[Omega]]Sqrt[(a^2+b^2)^(1/4)-Sqrt[(Sqrt[a^2+b^2]+a)/2.]]}/Sqrt[2.],
Abs[p^2-pc2]<=10^-15,
{Sqrt[Abs@b](Sqrt[Sqrt[2.]+1]+I Sign[\[Omega]]Sqrt[Sqrt[2.]-1])
}/Sqrt[2.],
p^2>pc2,
{Sqrt[Sqrt[a^2+b^2]-a]+I Sign[\[Omega]]Sqrt[Sqrt[a^2+b^2]+a],
Sqrt[(a^2+b^2)^(1/4)+Sqrt[(Sqrt[a^2+b^2]-a)/2.]]+I Sign[\[Omega]]Sqrt[(a^2+b^2)^(1/4)-Sqrt[(Sqrt[a^2+b^2]-a)/2.]]}/Sqrt[2.]
]
];

{q\[Alpha],q\[Alpha]s}=Q[#[[1]],#[[2]]]&/@{q\[Alpha]2,d/c[11]}


{q\[Beta],q\[Beta]s}


d2=Sqrt[(\[Sigma]0+\[Delta])/(\[Sigma]0+S\[Alpha])];
d3=2\[Beta]0^2 (\[Sigma]0+1/2 (S\[Alpha]+\[Delta]))/(\[Sigma]0+\[Delta]);
d4=Sqrt[(\[Sigma]0-p^2 \[Beta]0^2 (\[Sigma]0+S\[Beta]))/((1-p^2 \[Beta]0^2 (1+S\[Beta]))(\[Sigma]0+\[Delta]))];
d5=(\[Sigma]0-2p^2 \[Beta]0^2 (\[Sigma]0+1/2 (S\[Beta]+\[Delta])))/(\[Sigma]0+\[Delta]);
d1=1/Sqrt[p^2 d3+d5];
ds=Sqrt[d];

L1=d1{
{d2 q\[Alpha]s/ds,1/d4 p/ds/q\[Beta]s},
{d3 d2 p ds q\[Alpha]s,-d5/d4 ds/q\[Beta]s}
};

L2=d1{
{d5/d2 ds/q\[Alpha]s,d3 d4 p ds q\[Beta]s},
{1/d2 p/ds/q\[Alpha]s,-d4 q\[Beta]s/ds}
};
Developer`ToPackedArray/@{{q\[Alpha],q\[Beta]},L1,L2}
]
];
tlQLmatricesRadicals[index_Integer,p_]:=(Message[tlQLmatricesRadicals::nodef,index];
$Failed);
tlQLmatricesRadicals[index_,p_]:=(Message[tlQLmatricesRadicals::noint,index];
$Failed);
*)


(* ::Input::Initialization:: *)
SetAttributes[tlLayerReflectionTransmission,Listable];tlLayerReflectionTransmission[index1_Integer,index2_Integer,p_]/;(ListQ[tlLayer[index1]]&&ListQ[tlLayer[index2]]):=
Module[{L1t, L2t, L1b, L2b, C, D, CpD},
{L1t, L2t, L1b, L2b}=tlQLmatrices[index1,p][[2;;]]~Join~tlQLmatrices[index2,p][[2;;]];
C=Transpose[L2t].L1b;
D=Transpose[L1t].L2b;
CpD=Inverse[C+D];
{2CpD,(C-D).CpD}
];
tlLayerReflectionTransmission[index_Integer,p_]:=tlLayerReflectionTransmission[index, index+1,p];
tlLayerReflectionTransmission[index_Integer,index2_Integer,p_]:=(Message[tlLayerReflectionTransmission::nodef,index];
$Failed);
tlLayerReflectionTransmission[index_,index_,p_]:=(Message[tlLayerReflectionTransmission::noint,index];
$Failed);


(* ::Input::Initialization:: *)
tlReflectivity[indexSet:{_Integer..},thicknessSet:{_?NonNegative..}, p_, \[Omega]_]/;(Length@thicknessSet==Length@indexSet):=
Module[{f, foldList, CDmats, TRmats, qL,Emats},
qL=tlQLmatrices[indexSet,p,\[Omega]];
CDmats=Transpose@MapThread[{Transpose[#1[[2]]].(#2[[1]]),Transpose[#1[[1]]].(#2[[2]])}&,{qL[[;;,{2,3}]],RotateLeft@qL[[;;,{2,3}]]}];
TRmats=Transpose@MapThread[Block[{temp=Inverse[#1+#2]},{2temp,(#1-#2).temp}]&,CDmats];

Emats=DiagonalMatrix/@Exp[I \[Omega] thicknessSet (qL[[;;,1]])];

foldList=Reverse@Most@Transpose@{Sequence@@TRmats,Emats};

f=(#4.(#3+Transpose[#2].#1.Inverse[IdentityMatrix[2]+#3.#1].#2).#4)&;
Fold[f[#1,Sequence@@#2]&,{{0,0},{0,0}},foldList]//N
];


(* ::Input::Initialization:: *)
tlResponse[indexSet:{_Integer..},thicknessSet:{_?NonNegative..}, p_, \[Omega]_,zs_Real,zb_Real]/;(Length@thicknessSet==Length@indexSet):=
Module[{ref,qtop,L1top,L2top,Linv,vforce,\[CapitalSigma],S,uvec,bvec,StressDisplacement},
vforce={0,0,1/\[Omega],0};
{qtop,L1top,L2top}=tlQLmatrices[indexSet[[1]],p,\[Omega]];
ref=Block[{e=DiagonalMatrix[Exp[I \[Omega] qtop zb]]},e.tlReflectivity[indexSet,thicknessSet, p, \[Omega]].e];
Linv=-1/Sqrt[2] Transpose@Join[(I L2top)~Join~(-L1top),(I L2top)~Join~L1top,2];
\[CapitalSigma]=Linv.vforce{1,1,-1,-1};
S=DiagonalMatrix[Exp[I \[Omega] qtop~Join~(-qtop) zs]].\[CapitalSigma];
uvec=ref.S[[3;;4]]-S[[1;;2]];
bvec=1/Sqrt[2] {Sequence@@(I L1top.uvec),Sequence@@(L2top.uvec)};
StressDisplacement=bvec \[Omega]^2 {1/\[Omega],-1,1,1/\[Omega]}//N;
StressDisplacement[[{1,4}]]
];


(* ::Input::Initialization:: *)
tlResponseVforce[indexSet:{_Integer..},thicknessSet:{_?NonNegative..}, p_, \[Omega]_,zs_Real,zb_Real]/;(Length@thicknessSet==Length@indexSet):=
Module[{ref,qtop,L1top,L2top,Linv,vforce,\[CapitalSigma],S1,S2,uvec,bvec,StressDisplacement,e},
(*VERTICAL FORCE ONLY. FOR GENERAL SOURCE, USE tlResponse*)
vforce={1/\[Omega],0};
{qtop,L1top,L2top}=tlQLmatrices[indexSet[[1]],p,\[Omega]];
e[z_]:=N[DiagonalMatrix[Exp[I \[Omega] qtop z]]];
ref=e[zb].tlReflectivity[indexSet,thicknessSet, p, \[Omega]].e[zb];
\[CapitalSigma]=1/Sqrt[2]Transpose[L1top].vforce;
uvec=ref.(e[-zs].\[CapitalSigma]);
(*
<<WITH THE SOURCE TERM>>
\[CapitalSigma]=1/Sqrt[2]Transpose[L1top].vforce;
{S1,S2}={e[zs].\[CapitalSigma],e[-zs].\[CapitalSigma]};
uvec=ref.S2-S1;
*)
bvec=1/Sqrt[2] {Sequence@@(I L1top.uvec),Sequence@@(L2top.uvec)};
StressDisplacement=bvec \[Omega]^2 {1/\[Omega],-1,1,1/\[Omega]}//N;
(*p StressDisplacement[[{1,4}]]{BesselJ[0,p \[Omega] r],BesselJ[1,p \[Omega] r]} ADD r to the input*)
StressDisplacement[[{1,4}]]
];


(* ::Input::Initialization:: *)
tlDHTR[u_?VectorQ,\[Nu]_,r_,rmax_,BesselVector_?VectorQ]/;(Length@u==Length@BesselVector):=Module[{dhtlist,scorrection},
scorrection=Sin[#]/#&/@(BesselVector Pi/Last[BesselVector]);
dhtlist=2/(rmax^2 N[BesselJ[\[Nu]+1,#]^2])&/@BesselVector;
(scorrection dhtlist u).(N[BesselJ[\[Nu],# r/rmax]]&/@BesselVector)
];


(* ::Input::Initialization:: *)
tlPlotSeismic[label_?StringQ,u_,\[Omega]_,offsets_,ampl_,\[Omega]I_,model_,dlist_,{zs_,zb_},output_:0,fRicker_:25,tmax_,code_:0]/;(Length@model==Length@dlist):=
Module[{rickerSpectrum,dispfield,reflSource,reflTX,t,mm,HW,pcrit},
rickerSpectrum[w_,wp_]:=N[(2 w^2)/(Sqrt[\[Pi]] wp^3) Exp[-(w^2/wp^2)]];

dispfield=u;
HW=HannWindow[\[Omega]/2/Max@\[Omega]];
reflSource=HW dispfield(rickerSpectrum[#,2Pi fRicker]&/@\[Omega]);
(*
t=With[{len=Length@\[Omega],df=Abs[(\[Omega][[2]]-\[Omega][[1]])/(2Pi)]},Range[0,len]/(df (len))];
reflTX=Table[
Re[Exp[\[Omega]I t]InverseFourier[{0+0. I}~Join~(HW reflSource[[;;,i]])]],{i,1,Dimensions[reflSource][[2]]}];*)

t=With[{len=Length@\[Omega],df=(\[Omega][[2]]-\[Omega][[1]])/(2Pi)},Range[0,2len]/(df (2len))];
reflTX=Table[
Exp[\[Omega]I t]InverseFourier[{0+0. I}~Join~(reflSource[[;;,i]])~Join~(Reverse@Conjugate@reflSource[[;;,i]])],{i,1,Dimensions[reflSource][[2]]}];

If[Head@output==Symbol,Print["{data, t} written to "<>ToString[output]];output={reflTX,t};];

reflTX=ampl Reverse[Re@reflTX/Max@Abs@reflTX,2]\[Transpose];

mm=MinMax@reflTX;
pcrit=Min@Re@Flatten[({#[[1]],#[[4]]}/#[[-1]])^(-1/2)&/@tlGetStiffness[DeleteDuplicates@model]];
(*Print[mm];*)
Module[{pr={0,tmax},TToptns,is=500,tt=tlTravelTimes[model,dlist,#,zs,zb,code]&,lbls={{"t, s",None},{None,"x, km"}},cfg=Blend[{Black,White},(#-mm[[1]])/2]&,cfrb=Blend[{Blue,White,Red},(#+1)/2]&,colors={Cyan,Black,Green,Orange,Magenta}},
TToptns={PlotLabel->label,AspectRatio->1,Frame->True,FrameTicks->{{Range[0,Max@t,.5],None},{None,Range[0,Max@offsets,.25]}},ScalingFunctions->{Identity,"Reverse"},PlotRange->{MinMax@offsets,-pr},ImageSize->is,FrameTicksStyle->Opacity[0],FrameLabel->lbls,LabelStyle->Directive[Opacity[0],FontSize->20,Bold,Black],FrameStyle->Opacity[0]};

Overlay[{
ArrayPlot[reflTX,AspectRatio->1,Frame->True,FrameLabel->RotateRight@lbls,LabelStyle->Directive[FontSize->20,Bold,Black],DataReversed->{True,False},DataRange->{MinMax@(offsets),MinMax@t},FrameTicks->{{Range[0,Max@t,.5],None},{None,Range[0,Max@offsets,.25]}},PlotRange->{MinMax@(offsets),Max@t-Reverse@pr,Full},ImageSize->is,ColorFunction->cfrb,ColorFunctionScaling->False,PlotLabel->label,
PlotLegends->If[code==0,Placed[LineLegend[Directive[Thick,#]&/@colors,{"P, S","PP","SS","PS","SP"},LegendLayout->"Row"],Below],Placed[LineLegend[{Directive[Thick,Dashed]},{code},LegendLayout->"Row"],Below]]



],
ParametricPlot[Evaluate@tt[p][[2]],{p,0,1},PlotStyle->Directive[Dashed,Opacity[0],First@colors],#]&@TToptns,

If[code==0,
Sequence@@Table[
ParametricPlot[Evaluate@tt[p][[1]][[i]],{p,0,.99pcrit},PlotStyle->Directive[Dashed,Thick,colors[[i+1]]],#]&@TToptns,{i,1,4}
],
ParametricPlot[Evaluate@tt[p],{p,0,.99pcrit},PlotStyle->Directive[Dashed,Thick,Black],#]&@TToptns
]
}]

]
];


(* ::Input::Initialization:: *)
tlPlotDiff[label_?StringQ,delta_,offsets_,t_,ampl_,pr_,gain_,as_:1]:=Module[{reflTX,mm},
reflTX=({t[[2]]}~Join~Rest@t)^gain ampl Reverse[ Re@delta/Max@Abs@delta,2]\[Transpose];
mm=MinMax@reflTX;
Module[{TToptns,is=500,tt=tlTravelTimes[layers,depth,#,.4,.8]&,lbls={{"t, s",None},{None,"x, km"}},cfg=Blend[{Black,White},(#-mm[[1]])/2]&,cfrb=Blend[{Blue,White,Red},(#+1)/2]&,colors={Cyan,Black,Green,Orange,Magenta}},
TToptns={PlotLabel->label,AspectRatio->1,Frame->True,FrameTicks->{{Range[0,Max@t,.5],None},{None,Range[0,Max@offsets,.25]}},ScalingFunctions->{Identity,"Reverse"},PlotRange->{MinMax@offsets,-pr},ImageSize->is,FrameTicksStyle->Opacity[0],FrameLabel->lbls,LabelStyle->Directive[Opacity[0],FontSize->15,Bold,Black],FrameStyle->Opacity[0]};

ArrayPlot[reflTX,PlotLabel->label,AspectRatio->as,Frame->True,FrameLabel->RotateRight@lbls,LabelStyle->Directive[FontSize->20,Bold,Black],DataReversed->{True,False},DataRange->{MinMax@(offsets),MinMax@t},FrameTicks->{{Range[0,Max@t,.5],None},{None,Range[0,Max@offsets,.25]}},PlotRange->{MinMax@(offsets),Max@t-Reverse@pr,Full},ImageSize->is,ColorFunction->cfrb,ColorFunctionScaling->False]

]

]


(* ::Input::Initialization:: *)
tlPlotTrace[u_,\[Omega]_,\[Omega]I_,clr_:Black,fRicker_:25,pr_,is_,label_,prv_:Full]:=
Module[{rickerSpectrum,dispfield,reflSource,reflTX,t,mm},
rickerSpectrum[w_,wp_]:=(2 w^2)/(Sqrt[\[Pi]] wp^3) Exp[-(w^2/wp^2)];

dispfield=u;
t=With[{len=Length@\[Omega],df=(\[Omega][[2]]-\[Omega][[1]])/(2Pi)},Range[0,2len]/(df (2len))];

reflSource=dispfield(rickerSpectrum[#,2Pi fRicker]&/@\[Omega]);
reflTX=
Exp[\[Omega]I t]InverseFourier[{0+0. I}~Join~(reflSource)~Join~(Reverse@Conjugate@reflSource)];
Module[{lbls={{"Amplitude",None},{"t, s",label}},optns},

optns={LabelStyle->Directive[Black,Bold,FontSize->30],PlotStyle->Directive[Thick,clr],PlotRange->{pr,prv},ImageSize->is,Frame->True,AspectRatio->1/2,FrameLabel->lbls,GridLines->Automatic};

ListLinePlot[{t,reflTX}\[Transpose],#]&@optns
]
];


(* ::Input::Initialization:: *)
tlTravelTimes[indexSet:{_Integer..},thicknessSet:{_?NonNegative..}, p_, zs_Real,zb_Real,code_:0]/;(Length@thicknessSet==Length@indexSet):=
Module[{dqp2,dqs2,\[Alpha]0,\[Beta]0,\[Eta],\[Delta],\[Sigma]0,dq,xPS,tUD,Q,dQ,c11,c13,c33,c44,c55,d,xDct,tDct,PP,SS,PS,SP,Refl},
{c11, c13, c33, c44, d} = Re@Transpose[tlGetStiffness/@indexSet];
\[Alpha]0=Sqrt[c33/d];
\[Beta]0=Sqrt[c44/d];
\[Sigma]0=1-c44/c33;
\[Delta]=(c13-c33+2c44)/c33;
\[Eta]=(c11 c33-(c13+2c44)^2)/(2c33^2);
 
dq=Sqrt[\[Sigma]0^2+4 p^4 \[Alpha]0^2 \[Eta] (\[Alpha]0^2 \[Eta]+2 \[Beta]0^2 (\[Delta]+\[Sigma]0))+4 p^2 (-\[Alpha]0^2 \[Eta] \[Sigma]0+\[Beta]0^2 \[Delta] (\[Delta]+\[Sigma]0))];
dqp2=Re[-((2 p (\[Alpha]0^2 \[Eta] (dq+2 p^2 \[Alpha]0^2 \[Eta]-\[Sigma]0)+\[Beta]0^2 (dq+dq \[Delta]+(\[Delta]+4 p^2 \[Alpha]0^2 \[Eta]) (\[Delta]+\[Sigma]0))))/(\[Beta]0^2 dq))];
dqs2=Re[p (-4-4 \[Delta]-(4 \[Alpha]0^2 \[Eta])/\[Beta]0^2)-dqp2];
 
Q=Re[tlQLmatrices[indexSet,p][[;;,1]]];
 
dQ=Re[Transpose[{dqp2,dqs2}]/2/Q];
Block[{z={zb}~Join~Rest@thicknessSet},
xPS=-dQ z;tUD=xPS p +Q z;];

(*DIRECT*)
xDct=-dQ[[1]]zs;
tDct=xDct p+Q[[1]]zs;

 
Block[{dir=Transpose@{xDct,tDct},seq,layers},
If[code==0,

{PP,SS,PS}=Accumulate/@Transpose/@{
{2xPS[[;;,1]],2tUD[[;;,1]]},
{2xPS[[;;,2]],2tUD[[;;,2]]},
{Total[xPS,{2}],Total[tUD,{2}]}
};
 
{Transpose[(-{Sequence@@dir,Sequence@@dir}+#)&/@Transpose[{PP,SS,PS,PS}]],dir}

 ,

If[StringLength@code==1,
Switch[code,"P",dir[[1]],"S",dir[[2]]],

layers=ToExpression@DeleteCases[StringSplit[code,_?LetterQ],""];
seq=ToExpression@StringReplace[#,{"P"->"1","S"->"2"}]&/@StringSplit[#,""]&/@StringSplit[code,_?DigitQ];
Total[Flatten[Table[{xPS[[layers[[i]],#]],tUD[[layers[[i]],#]]}&/@(seq[[i]]),{i,Length@layers}],1],{1}]-dir[[First@First@seq]]
	

]]]
];


(* ::Input::Initialization:: *)
Options[tlVisualiseModel]={ColorSchemeNumber->3, SamplingRate->1000};
tlVisualiseModel[indexSet:{_Integer..},thicknessSet:{_?Positive..}, nTimeSamples_Integer, opts:OptionsPattern[]]/;(Length@thicknessSet==Length@indexSet):=Module[{layers, colourRules,dt, sr},
layers=Sort@DeleteDuplicates@indexSet;
If[IntegerQ@OptionValue[ColorSchemeNumber],
colourRules=Thread[layers->ColorData[OptionValue[ColorSchemeNumber]]/@Range@Length@layers],
colourRules=Thread[layers->(ColorData[3]/@Range@Length@layers)]
];
If[IntegerQ@OptionValue[SamplingRate]&&OptionValue[SamplingRate]>499,
sr=OptionValue[SamplingRate],
sr=1000
];
dt=Block[{c, d},
{c[33],c[44],d}=tlLayer[#1][[3;;]];
{#2/Sqrt[c[33]/d],#2/Sqrt[c[44]/d]}
]&;
ArrayPlot[Transpose[ConstantArray[Join@@MapThread[ConstantArray[#1,Round[sr dt[#1,#2][[1]],1]]&,Reverse/@{Most@indexSet, Most@thicknessSet}],5]],ColorRules->colourRules]
];


(* ::Input::Initialization:: *)
Options[outerFunction] = 
   {option->value}; 
outerFunction[args__, opts:OptionsPattern[]] := privatefunction[args, 
OptionValue[option]
];


(* ::Input::Initialization:: *)
End[]


(* ::Input::Initialization:: *)
EndPackage[]
