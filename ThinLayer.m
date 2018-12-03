(* ::Package:: *)

(* ::Input::Initialization:: *)
BeginPackage["ThinLayer`"]


(* ::Input::Initialization:: *)
(* Usage Messaging *)
tlSetLayer::usage="tlSetLayer[layerIndex, {c11, c13, c33, c44, \[Rho]}] creates a 'tlLayer[layerIndex]' object instance with assigned numerical values for cijs.";
tlGetStiffness::usage="tlGetStiffness[layerIndex] outputs a vector of cijs for a given layer.";
tlCijFromThomsen::usage="tlCijFromThomsen[vp0,vs0,\[Epsilon],\[Delta],\[Rho]] outputs {c11, c13, c33, c44, \[Rho]} for a set of Thomsen parameters.";
tlVTISquirtModuli::usage="tlVTISquirtModuli[\[Epsilon],\[Epsilon]f,\[Alpha],\[Tau]0, \[Omega], lenrat][\[Lambda],\[Mu], \[Phi], Kf] returns {c11, c13, c33, c44} for a given frequency.";
tlQLmatrices::usage="tlQmatices[i, p] returns the numbers qa, qb and linear transformations L1, L2 decomposing the wave field of horizontal slowness p, in a layer with index i. The output is in the form {{qa, qb}, L1, L2}";
tlLayerReflectionTransmission::usage="tlLayerReflectionTransmission[i, j, p] calculates transmission T and reflection R matrices for an interface between pre-defined layers i and j for horizontal slowness p. It can also be called as tlLayerReflectionTransmission[i, p] in which case an interface between layers i, i+1 is assumed. The output is of the form {T, R}."; 
tlReflectivity::usage="tlReflectivity[{i,j, ...},{di, dj, ...}, slowness, angularFreq] calculates the reflection response matrix RR of a sequence of pre-initiated layers labelled {i,j,...} with thicknesses {di, dj, ...} where the first and last layers are halfspaces and the first and last sequence is ignored.";
tlResponse::usage="tlResponse[{i,j, ...},{di, dj, ...}, slowness, angularFreq,zs,zb] calculates the stress-displacement vector for a stack of layers bounded by two half-spaces. source and receiver are located in the top half-space: source depth: zs, receiver depth: 0, first boundary: zb. the vetor returned is in frequency-slowness domain: {Uz, Srz, Szz, Ur}";
tlResponseVforce::usage="tlResponseVforce[{i,j, ...},{di, dj, ...}, slowness, angularFreq,zs,zb] calculates the stress-displacement vector for a stack of layers bounded by two half-spaces. source and receiver are located in the top half-space: source depth: zs, receiver depth: 0, first boundary: zb. the vetor returned is in frequency-slowness domain: {Uz, Srz, Szz, Ur}. Source type - vertical force";
tlDiscreteHankelTransform::usage="tlDiscreteHankelTransform[{f(k1), f(k2), f(k3), ..., f(kmax)}, ord, kmax] returns {xi, H(fi)} as pairs of offsets xi and the Hankel transform H(fi) of order ord of a function f evaluated at k1, ... k(max)";
tlDHT2::usage="tlDHT2[N, f[k], ord, kmax] returns {xi, H(f[k])i} as pairs of offsets xi and the Hankel transform H(f[k]) of order ord of a function f[k] passed in as a function";
tlDHTp::usage="tlDHT2[N, f[k], ord, kmax, \[Omega]] returns {xi, H(f[k])i} as pairs of offsets xi and the Hankel transform H(f[k]) of order ord of a function f[k] passed in as a function. \[Omega] is passed when the transform is applied in slowness domain.";
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
\[Phi]c=(4 \[Pi] \[Epsilon] \[Alpha])/3;
\[Phi]f=(4 \[Pi] \[Epsilon]f \[Alpha])/3;
\[Phi]p=\[Phi]-\[Phi]c-\[Phi]f;
Kp=(4\[Mu])/(3 Kf);
Kc=\[Sigma]c/Kf;
\[Nu]=\[Lambda]/(2(\[Lambda]+\[Mu]));(* Poisson's ratio *)
\[Sigma]c=(\[Pi] \[Mu] \[Alpha])/(2(1-\[Nu])); (* Crack stiffness *)

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
Module[{c11,c33,c13,c44},
c11=(1+2\[Epsilon])vp0^2 \[Rho];
c33=vp0^2 \[Rho];
c13=(-vs0^2+Sqrt[(vp0-vs0) (vp0+vs0) (-vs0^2+vp0^2 (1+2 \[Delta]))]) \[Rho];
c44=vs0^2 \[Rho];
{c11,c13,c33,c44,\[Rho]}//N
]



(* ::Input::Initialization:: *)
SetAttributes[tlQLmatrices,Listable];
tlQLmatrices[index_Integer,p_]/;ListQ[tlLayer[index]]:=
Block[{c, d},
{c[11],c[13],c[33],c[44],d}=tlLayer[index];
Module[
{L1,L2,d1,d2,d3,d4,d5,q\[Alpha],q\[Beta],\[Alpha]0,\[Beta]0,\[Eta],\[Delta],\[Sigma]0,S\[Alpha],S\[Beta],R,R1,R2,req1,req2,imq1,imq2},
\[Alpha]0=Sqrt[c[33]/d];
\[Beta]0=Sqrt[c[44]/d];
\[Sigma]0=1-c[44]/c[33];
\[Delta]=(c[13]-c[33]+2c[44])/c[33];
\[Eta]=(c[11] c[33]-(c[13]+2c[44])^2)/(2c[33]^2);
R1=2(1-p^2 \[Beta]0^2)(\[Delta]+2p^2 \[Alpha]0^2 \[Eta])^2;
R2=\[Sigma]0+2p^2 \[Beta]0^2 \[Delta]-2p^2 \[Alpha]0^2 (1-2p^2 \[Beta]0^2)\[Eta];
R=R1/(R2+Sqrt[R2^2+2p^2 \[Beta]0^2 R1]);
S\[Alpha]=2\[Delta]+2p^2 \[Alpha]0^2 \[Eta]+R//Chop;
S\[Beta]=2(1-p^2 \[Beta]0^2) \[Alpha]0^2/\[Beta]0^2 \[Eta]-R//Chop;
req1=Re[1/\[Alpha]0^2-p^2-p^2 S\[Alpha]];
req2=Re[1/\[Beta]0^2-p^2-p^2 S\[Beta]];
imq1=Im[1/\[Alpha]0^2-p^2-p^2 S\[Alpha]];
imq2=Im[1/\[Beta]0^2-p^2-p^2 S\[Beta]];
imq1=Sign[imq1]imq1;
imq2=Sign[imq2]imq2;
q\[Alpha]=Sqrt[req1+I imq1]//Chop;
q\[Beta]=Sqrt[req2+I imq2]//Chop;
d2=Sqrt[(\[Sigma]0+\[Delta])/(\[Sigma]0+S\[Alpha])];
d3=2\[Beta]0^2 (\[Sigma]0+1/2 (S\[Alpha]+\[Delta]))/(\[Sigma]0+\[Delta]);
d4=Sqrt[(\[Sigma]0-p^2 \[Beta]0^2 (\[Sigma]0+S\[Beta]))/((1-p^2 \[Beta]0^2 (1+S\[Beta]))(\[Sigma]0+\[Delta]))];
d5=(\[Sigma]0-2p^2 \[Beta]0^2 (\[Sigma]0+1/2 (S\[Beta]+\[Delta])))/(\[Sigma]0+\[Delta]);
d1=1/Sqrt[p^2 d3+d5];
L1=d1{
{d2 Sqrt[q\[Alpha]/d],1/d4 p/Sqrt[d q\[Beta]]},
{d3 d2 p Sqrt[d q\[Alpha]],-d5/d4 Sqrt[d/q\[Beta]]}
}//Chop;
L2=d1{
{d5/d2 Sqrt[d/q\[Alpha]],d3 d4 p Sqrt[d q\[Beta]]},
{1/d2 p/Sqrt[d q\[Alpha]],-d4 Sqrt[q\[Beta]/d]}
}//Chop;
Developer`ToPackedArray/@{{q\[Alpha],q\[Beta]},L1,L2}
]
];
tlQLmatrices[index_Integer,p_]:=(Message[tlQLmatrices::nodef,index];
$Failed);
tlQLmatrices[index_,p_]:=(Message[tlQLmatrices::noint,index];
$Failed);


(* ::Input::Initialization:: *)
SetAttributes[tlLayerReflectionTransmission,Listable];tlLayerReflectionTransmission[index1_Integer,index2_Integer,p_]/;(ListQ[tlLayer[index1]]&&ListQ[tlLayer[index2]]):=
Module[{L1t, L2t, L1b, L2b, C, D, CpD},
{L1t, L2t, L1b, L2b}=tlQLmatrices[index1,p][[2;;]]~Join~tlQLmatrices[index2,p][[2;;]];
C=Transpose[L2b].L1t;
D=Transpose[L1b].L2t;
CpD=Transpose@Inverse[C+D];
{2CpD, Transpose[C-D].CpD}
];
tlLayerReflectionTransmission[index_Integer,p_]:=tlLayerReflectionTransmission[index, index+1,p];
tlLayerReflectionTransmission[index_Integer,index2_Integer,p_]:=(Message[tlLayerReflectionTransmission::nodef,index];
$Failed);
tlLayerReflectionTransmission[index_,index_,p_]:=(Message[tlLayerReflectionTransmission::noint,index];
$Failed);


(* ::Input::Initialization:: *)
tlReflectivity[indexSet:{_Integer..},thicknessSet:{_?NonNegative..}, p_, \[Omega]_]/;(Length@thicknessSet==Length@indexSet):=
Module[{f, foldList, CDmats, TRmats, qL},
qL=tlQLmatrices[indexSet,p];
CDmats=Transpose@MapThread[
{Transpose[#1[[2]]].#2[[1]],Transpose[#1[[1]]].#2[[2]]}&,
{RotateLeft@qL[[;;,{2,3}]],qL[[;;,{2,3}]]}];
TRmats=Transpose@MapThread[Block[{temp=Transpose@Inverse[#1+#2]},{2temp,Transpose[(#1-#2)].temp}]&,CDmats];
foldList=Reverse@Most@Transpose@{Sequence@@TRmats,DiagonalMatrix/@Exp[I \[Omega] thicknessSet({{0.,0.}}~Join~Rest@Most@qL[[;;,1]]~Join~{{0.,0.}})]};
f=(#4.(#3+Transpose[#2].#1.Inverse[IdentityMatrix[2]+#3.#1].#2).#4)&;
Fold[f[#1,Sequence@@#2]&,{{0,0},{0,0}},foldList]//N//Chop
];


(* ::Input::Initialization:: *)
tlResponse[indexSet:{_Integer..},thicknessSet:{_?NonNegative..}, p_, \[Omega]_,zs_Real,zb_Real]/;(Length@thicknessSet==Length@indexSet):=
Module[{ref,qtop,L1top,L2top,Linv,vforce,\[CapitalSigma],S,uvec,bvec,StressDisplacement},
vforce={0,0,-1/\[Omega],0};
{qtop,L1top,L2top}=tlQLmatrices[indexSet[[1]],p];
ref=Block[{e=DiagonalMatrix[Exp[I \[Omega] qtop zb]]},e.tlReflectivity[indexSet,thicknessSet, p, \[Omega]].e];
Linv=-1/Sqrt[2] Transpose@Join[(I L2top)~Join~(-L1top),(I L2top)~Join~L1top,2];
\[CapitalSigma]=Linv.vforce{1,1,-1,-1};
S=DiagonalMatrix[Exp[I \[Omega] qtop~Join~(-qtop) zs]].\[CapitalSigma];
uvec=ref.S[[3;;4]]-S[[1;;2]];
bvec=1/Sqrt[2] {Sequence@@(I L1top.uvec),Sequence@@(-L2top.uvec)};
StressDisplacement=bvec \[Omega]^2 {1/\[Omega],-1,1,1/\[Omega]}//N;
StressDisplacement[[{1,4}]]
];


(* ::Input::Initialization:: *)
tlResponseVforce[indexSet:{_Integer..},thicknessSet:{_?NonNegative..}, p_, \[Omega]_,zs_Real,zb_Real,r_Real]/;(Length@thicknessSet==Length@indexSet):=
Module[{ref,qtop,L1top,L2top,Linv,vforce,\[CapitalSigma],S1,S2,uvec,bvec,StressDisplacement,e},
e[z_]:=N[DiagonalMatrix[Exp[I \[Omega] qtop z]]];
(*VERTICAL FORCE ONLY. FOR GENERAL SOURCE, USE tlResponse*)
vforce={-1/\[Omega],0};
{qtop,L1top,L2top}=tlQLmatrices[indexSet[[1]],p];
ref=e[zb].tlReflectivity[indexSet,thicknessSet, p, \[Omega]].e[zb];
\[CapitalSigma]=1/Sqrt[2]Transpose[L1top].vforce;
{S1,S2}={e[zs].\[CapitalSigma],e[-zs].\[CapitalSigma]};
uvec=ref.S2-S1;
bvec=1/Sqrt[2] {Sequence@@(I L1top.uvec),Sequence@@(-L2top.uvec)};
StressDisplacement=bvec \[Omega]^2 {1/\[Omega],-1,1,1/\[Omega]}//N;
(*StressDisplacement[[{1,4}]];*)
p StressDisplacement[[{1,4}]]{BesselJ[0,p \[Omega] r],BesselJ[1,p \[Omega] r]}
];


(* ::Input::Initialization:: *)
tlDiscreteHankelTransform[f1_?VectorQ,p_:0.,rmax_]:=Module[{Np,a,aNp1,rv,uv,res,umax,T,J,F1,F2},
Np=Length@f1;
a=Developer`ToPackedArray@Table[N[BesselJZero[p,n]],{n,1,Np}];
aNp1=N[BesselJZero[p,Np+1]];
umax=aNp1/(2. Pi rmax);
rv=a/(2. Pi umax);
uv=a/(2. Pi rmax);
J=Developer`ToPackedArray@Abs@Table[BesselJ[p+1,s],{s,a}];
T=Developer`ToPackedArray@(Table[BesselJ[p,s/(2. Pi rmax umax)],{s,TensorProduct[a,a]}]/(TensorProduct[J,J] Pi rmax umax));
F1=(f1 rmax)/J;
F2=Developer`ToPackedArray[T.F1];
Transpose[{uv+0.I,J/umax F2}]
];


(* ::Input::Initialization:: *)
tlDHT2[Np_?IntegerQ,ReflFun_,p_:0.,kmax_]:=Module[{a,aNp1,kvec,rvec,rmax,T,J,F1,F2,f1,f2,pi},
pi=1. ;(* 1 or 2. Pi *)
a=Developer`ToPackedArray@Table[N[BesselJZero[p,n]],{n,1,Np}];
aNp1=N[BesselJZero[p,Np+1]];
kvec=a kmax/aNp1;
rvec=a/(pi kmax);
rmax=aNp1/(pi kmax);
J=Developer`ToPackedArray@Abs@Table[N[BesselJ[p+1,s]],{s,a}];
T=Developer`ToPackedArray@(Table[2N[BesselJ[p,s/(pi rmax kmax)]],{s,TensorProduct[a,a]}]/(TensorProduct[J,J]pi rmax kmax));
f2=ReflFun[#]&/@kvec;
F2=kmax f2/J;
F1=Developer`ToPackedArray[T.F2];
f1=F1 J/rmax;
Transpose[{rvec,f1}]
];


(* ::Input::Initialization:: *)
tlDHTp[Np_?IntegerQ,ReflFun_,p_:0.,pmax_,\[Omega]_:1.]:=Module[{a,aNp1,pvec,rvec,rmax,T,J,F1,F2,f1,f2,pi},
pi=1. ;(* 1 or 2. Pi *)
a=Developer`ToPackedArray@Table[N[BesselJZero[p,n]],{n,1,Np}];
aNp1=N[BesselJZero[p,Np+1]];
pvec=a pmax/aNp1;
rvec=1/\[Omega] a/(pi pmax);
rmax=1/\[Omega] aNp1/(pi pmax);
J=Developer`ToPackedArray@Abs@Table[N[BesselJ[p+1,s]],{s,a}];
T=Developer`ToPackedArray@(Table[2N[BesselJ[p,s/(pi rmax pmax)]],{s,TensorProduct[a,a]}]/(TensorProduct[J,J]pi rmax pmax));
f2=ReflFun[#]&/@pvec;
F2=pmax f2/J;
F1=Developer`ToPackedArray[T.F2];
f1=F1 J/rmax;
Transpose[{rvec,f1}]
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
