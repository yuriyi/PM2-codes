(* ::Package:: *)

(* ::Input::Initialization:: *)
BeginPackage["ThinLayer`"]


(* ::Input::Initialization:: *)
(* Usage Messaging *)
tlSetLayer::usage="tlSetLayer[index, {c11, c13, c33, c44, \[Rho]}] creates a 'tlLayer[index]' object instance with numerical values for cijs.";
tlQmatrices::usage="tlQmatices[i, p] calculates the Q matrices corresponding to a pre-defined layer with index i for a horizontal slowness p";
tlLmatrices::usage="tlLmatices[i, p] calculates the L matrices corresponding to a pre-defined layer with index i for a horizontal slowness p";
tlLayerReflectionTransmission::usage="tlQLmatices[i, j, p] calculates the reflection and transmission coefficients for an interface at layer \!\(\*
StyleBox[\"i\",\nFontSlant->\"Italic\"]\) and \!\(\*
StyleBox[\"j\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)for horizontal slowness \!\(\*
StyleBox[\"p\",\nFontSlant->\"Italic\"]\). Can also be called as tlQLmatices[i, p] in which case the coefficients at the interface between layers \!\(\*
StyleBox[\"i\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"i\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"+\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"1\",\nFontSlant->\"Italic\"]\) is calculated"; 
tlLayerPhase::usage="";
tlReflectivity::usage=""; 
(* Error Reporting *)
tlQmatrices::nodef=tlLayerReflectionTransmission::nodef="Layer `1` not defined. Define elastic constants for layer `1` using tlSetLayer";
tlQmatrices::noint=tlLayerReflectionTransmission::noint="Argument `1` refers to layer label. Integer expected";


(* ::Input::Initialization:: *)
Begin["Private`"]


(* ::Input::Initialization:: *)
tlSetLayer[index_Integer,{c11_?NumericQ, c13_?NumericQ, c33_?NumericQ, c44_?NumericQ, \[Rho]_Real}]:=
tlLayer[index] = {c[11]:>N@c11, c[13]:>N@c13,c[33]:>N@c33,c[44]:>N@c44, d:>N@\[Rho] }; 


(* ::Input::Initialization:: *)
SetAttributes[tlQmatrices,Listable];
tlQmatrices[index_Integer,p_]/;ListQ[tlLayer[index]]:=
Block[{c, d},
Set@@@tlLayer[index];
Module[
{q\[Alpha],q\[Beta],\[Alpha]0,\[Beta]0,\[Eta],\[Delta],\[Sigma]0,S\[Alpha],S\[Beta],R,R1,R2,req1,req2,imq1,imq2},
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
req1=Re[1/\[Alpha]0^2-p^2-p^2 S\[Alpha]];
req2=Re[1/\[Beta]0^2-p^2-p^2 S\[Beta]];
imq1=Im[1/\[Alpha]0^2-p^2-p^2 S\[Alpha]];
imq2=Im[1/\[Beta]0^2-p^2-p^2 S\[Beta]];
imq1=Sign[imq1]imq1;
imq2=Sign[imq2]imq2;
q\[Alpha]=Sqrt[req1+I imq1]//Chop;
q\[Beta]=Sqrt[req2+I imq2]//Chop;
{q\[Alpha],q\[Beta], S\[Alpha], S\[Beta]}
]
];
tlQmatrices[index_Integer,p_]:=(Message[tlQLmatrices::nodef,index];
$Failed);
tlQmatrices[index_,p_]:=(Message[tlQLmatrices::noint,index];
$Failed);


(* ::Input::Initialization:: *)
SetAttributes[tlLmatrices,Listable];
tlLmatrices[index_Integer,p_]/;ListQ[tlLayer[index]]:=
Block[{c, d},
Set@@@tlLayer[index];
Module[
{L1,L2,d1,d2,d3,d4,d5,q\[Alpha],q\[Beta],\[Alpha]0,\[Beta]0,\[Eta],\[Delta],\[Sigma]0,S\[Alpha],S\[Beta]},
\[Alpha]0=Sqrt[c[33]/d];
\[Beta]0=Sqrt[c[44]/d];
\[Sigma]0=1-c[44]/c[33];
\[Delta]=(c[13]-c[33]+2c[44])/c[33];
\[Eta]=(c[11] c[33]-(c[13]+2c[44])^2)/(2c[33]^2);
{q\[Alpha],q\[Beta],S\[Alpha],S\[Beta]}=tlQmatrices[index,p];

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
{L1,L2}
]
];
tlLmatrices[index_Integer,p_]:=(Message[tlQLmatrices::nodef,index];
$Failed);
tlLmatrices[index_,p_]:=(Message[tlQLmatrices::noint,index];
$Failed);


(* ::Input::Initialization:: *)
SetAttributes[tlLayerReflectionTransmission,Listable];tlLayerReflectionTransmission[index1_Integer,index2_Integer,p_]/;(ListQ[tlLayer[index1]]&&ListQ[tlLayer[index2]]):=
Module[{L1t, L2t, L1b, L2b, C, D, CpD},
{L1t, L2t, L1b, L2b}=tlLmatrices[index1,p]~Join~tlLmatrices[index2,p];
C=Transpose[L2b].L1t;
D=Transpose[L1b].L2t;
CpD=Transpose@Inverse[C+D];
{2CpD, Transpose[C-D].CpD}//Chop
];
tlLayerReflectionTransmission[index_Integer,p_]:=tlLayerReflectionTransmission[index, index+1,p];
tlLayerReflectionTransmission[index_Integer,index2_Integer,p_]:=(Message[tlLayerReflectionTransmission::nodef,index];
$Failed);
tlLayerReflectionTransmission[index_,index_,p_]:=(Message[tlLayerReflectionTransmission::noint,index];
$Failed);


(* ::Input::Initialization:: *)
tlLayerPhase[indexSet:{_Integer..},thicknessSet:{_?Positive..}, p_, \[Omega]_]/;(Length@thicknessSet==Length@indexSet):=
Exp[I \[Omega] thicknessSet tlQmatrices[indexSet,p][[;;,;;2]]]//Chop;


(* ::Input::Initialization:: *)
tlReflectivity[indexSet:{_Integer..},thicknessSet:{_?Positive..}, p_, \[Omega]_]/;((Length@thicknessSet==Length@indexSet)&&(And@@((#==0)&/@(indexSet-Range[indexSet[[1]],Length[indexSet]])))):=
Module[{f, foldList},
foldList=Transpose@{tlLayerReflectionTransmission[Most@indexSet, p],tlLayerPhase[Most@indexSet, thicknessSet,p,\[Omega]]};
f=(#4.(#3+Transpose[#2].#1.Inverse[IdentityMatrix[2]+#3.#1].#2).#4)&;
Fold[f[#1,#2[[1,1]],#2[[1,2]],DiagonalMatrix@#2[[2]]]&,{{0,0},{0,0}},foldList]//Chop
];
tlReflectivity[indexSet:{_Integer..},thicknessSet:{_?Positive..}, p_, \[Omega]_]/;(Length@thicknessSet==Length@indexSet):=
Module[{f, foldList},
foldList=Transpose@{tlLayerReflectionTransmission[Most@indexSet, Most@RotateLeft@indexSet, p],tlLayerPhase[Most@indexSet, thicknessSet,p,\[Omega]]};
f=(#4.(#3+Transpose[#2].#1.Inverse[IdentityMatrix[2]+#3.#1].#2).#4)&;
Fold[f[#1,#2[[1,1]],#2[[1,2]],DiagonalMatrix@#2[[2]]]&,{{0,0},{0,0}},foldList]//Chop
];
(*
tlLayerReflectionTransmission[index_Integer,p_]:=tlLayerReflectionTransmission[index, index+1,p];
tlLayerReflectionTransmission[index_Integer,index2_Integer,p_]:=(Message[tlLayerReflectionTransmission::nodef,index];
$Failed);
tlLayerReflectionTransmission[index_,index_,p_]:=(Message[tlLayerReflectionTransmission::noint,index];
$Failed);*)


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
