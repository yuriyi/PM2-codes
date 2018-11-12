(* ::Package:: *)

(* ::Input::Initialization:: *)
BeginPackage["ThinLayer`"]


(* ::Input::Initialization:: *)
(* Usage Messaging *)
tlSetLayer::usage="tlSetLayer[layerIndex, {c11, c13, c33, c44, \[Rho]}] creates a 'tlLayer[layerIndex]' object instance with assigned numerical values for cijs.";
tlQLmatrices::usage="tlQmatices[i, p] returns the numbers qa, qb and linear transformations L1, L2 decomposing the wave field of horizontal slowness p, in a layer with index i. The output is in the form {{qa, qb}, L1, L2}";
tlLayerReflectionTransmission::usage="tlLayerReflectionTransmission[i, j, p] calculates transmission T and reflection R matrices for an interface between pre-defined layers i and j for horizontal slowness p. It can also be called as tlLayerReflectionTransmission[i, p] in which case an interface between layers i, i+1 is assumed. The output is of the form {T, R}."; 
tlReflectivity::usage="tlReflectivity[{i,j, ...},{di, dj, ...}, slowness, angularFreq] calculates the reflection response matrix RR of a sequence of pre-initiated layers labelled {i,j,...} with thicknesses {di, dj, ...} where the first and last layers are halfspaces and the first and last sequence is ignored.";
tlDiscreteHankelTransform::usage="tlDiscreteHankelTransform[{f(k1), f(k2), f(k3), ..., f(kmax)}, ord, kmax] returns {xi, H(fi)} as pairs of offsets xi and the Hankel transform H(fi) of order ord of a function f evaluated at k1, ... k(max)";
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
tlReflectivity[indexSet:{_Integer..},thicknessSet:{_?Positive..}, p_, \[Omega]_]/;(Length@thicknessSet==Length@indexSet):=
Module[{f, foldList, CDmats, TRmats, qL, \[Phi]},
qL=tlQLmatrices[indexSet,p];
CDmats=Transpose@MapThread[
{Transpose[#1[[2]]].#2[[1]],Transpose[#1[[1]]].#2[[2]]}&,
{RotateLeft@qL[[;;,{2,3}]],qL[[;;,{2,3}]]}];
TRmats=Transpose@MapThread[Block[{temp=Transpose@Inverse[#1+#2]},{2temp,Transpose[(#1-#2)].temp}]&,CDmats];
foldList=Reverse@Most@Transpose@{Sequence@@TRmats,DiagonalMatrix/@Exp[I \[Omega] thicknessSet({{0.,0.}}~Join~Rest@qL[[;;,1]])]};
f=(#4.(#3+Transpose[#2].#1.Inverse[IdentityMatrix[2]+#3.#1].#2).#4)&;
Fold[f[#1,Sequence@@#2]&,{{0,0},{0,0}},foldList]
];


(* ::Input::Initialization:: *)
tlDiscreteHankelTransform[f1_?VectorQ,p_:0.,rmax_]:=Module[{Np,a,aNp1,rv,uv,res,umax,T,J,F1,F2},
Np=Length@f1;
a=Developer`ToPackedArray@Table[N[BesselJZero[p,n]],{n,1,Np}];
aNp1=BesselJZero[p,Np+1];
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
