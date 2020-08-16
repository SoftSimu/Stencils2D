(************************************************************************)
(*                                                                      *)
(*       This file belongs to the paper "Stencils with isotropic        *)
(*       discretisation error for differential operators"               *)
(*                                                                      *)
(*       File created by: Michael Patra <patra@lorentz.leidenuniv.nl    *)
(*                                                                      *)
(************************************************************************)
(*                                                                      *)
(* This file contains the stencils in Mathematica form. The definitions *)
(* in this file appear in the same sequence as they appear in the       *)
(* paper. The coefficients as well as entire stencils are given. Both   *)
(* specific stencils and general stencils are provided.                 *)
(*                                                                      *)
(* The naming scheme of the symbols is as follows:                      *)
(*                                                                      *)
(*               ====== "Anisotropic", "Isotropic" or "General"         *)
(* cLaplacian2D4hIso21p                                                 *)
(* ||        | | |  === 21 point stencil                                *)
(* ||        | | === "Iso": isotropic stencil, "Aniso": anisotropic     *)
(* ||        | == O(h^4) stencil                                        *)
(* ||        == two-dimensional                                         *)
(* |========= "Laplacian", "Bilaplacian" or "GradLap"                   *)
(* = "c" means that the symbol is not a stencil but its coefficients    *)
(*   "e" gives the prefactor for the discretisation error term          *)
(*   if empty, this is the stencil itself                               *)
(*                                                                      *)
(* In addition to the specific stencils (e.g. "...Iso21p"), the lists   *)
(* of conditions for the coefficients of general stencils are provided. *)
(* The symbol for the conditions that arbitrary stencils have to fulfil *)
(* has a name ending on "Anisotropic", for isotropic stencils it ends   *)
(* on "Isotropic".                                                      *)
(*                                                                      *)
(* The arrangements of the coefficients within the stencils are defined *)
(* by symbols ending on "General". Lists of coefficients can be applied *)
(* to it to yield stencils, i.e. "Laplacian2D4hIso21p" is identical to  *)
(* "Laplacian2D4hGeneral /. cLaplacian2D4hIso21p"                       *)
(*                                                                      *)
(* Examples:                                                            *)
(*   cBilaplacian3D2hIsotropic                                          *)
(*       Gives a list of conditions that all isotropic O(h^2) stencils  *)
(*       for the three-dimensional Bilaplacian have to fulfil.          *)
(*   cBilaplacian3D2hAnisotropic                                        *)
(*       Gives a list of conditions that all O(h^2) stencils for the    *)
(*       three-dimensional Bilaplacian have to fulfil without imposing  *)
(*       the restriction that the discretisation error has to be        *)
(*       isotropic.                                                     *)
(*   eLaplacian3D2hIso21p                                               *)
(*       The isotropic O(h^2) 21-point stencil for the Laplacian in     *)
(*       three dimensions has a discretisation error that is equal to   *)
(*       h^2 times eLaplacian3D2hIso21p times the Bilaplacian.          *)
(*   Bilaplacian2D4hIso37p[f,x0,y0,h]                                   *)
(*       Computes the O(h^4) approximation, using an 37-point isotropic *)
(*       stencil, of the Bilaplacian in two dimensions of the function  *)
(*       f[x,y] at the point x0,y0 where a grid size of h is used.      *)
(*   Simplify[  Laplacian2D4hGeneral /. cLaplacian2D4hIso21p            *)
(*            - Laplacian2D4hIso21p ]                                   *)
(*       Checks whether the coefficients for the O(h^4) stencil for the *)
(*       two-dimensional Laplacian agree with the stencil itself. The   *)
(*       result of the simplification process will, of course, be zero. *)
(*                                                                      *)
(* If this file is read into Mathematica, a series of tests are run on  *)
(* these definition. All relevant properties are checked. The tests are *)
(* as follows:                                                          *)
(*  a) It is verified that the coefficients agree with the stencil,     *)
(*     e.g. that Laplacian2D4hGeneral /. cLaplacian2D4hIso21p indeed is *)
(*     identical to Laplacian2D4hIso21p                                 *)
(*  b) Then all specific anisotropic stencils are checked by plugging   *)
(*     in a polynomial of sufficiently high degree. From the difference *)
(*     between this result and the exact analytical result, a power     *)
(*     series expansion in h is done. For an O(h^p) stencil, the power  *)
(*     series expansion of order p-1 has to vanish.                     *)
(*  c) The same is now done for the specific isotropic stencils.        *)
(*     However, from the difference we also substract the analytical    *)
(*     value of the discretisation error times its prefactor defined by *)
(*     the symbol with the "e" prefix (see "nomenclature" above).       *)
(*     Consequently the order of the power series expansion is now      *)
(*     higher.                                                          *)
(*  d) Finally the lists with the general conditions for arbitrary and  *)
(*     and isotropic stencils are checked in the same way. The needed   *)
(*     simplifications of the equations are time-consuming and might    *)
(*     take several hours on older workstations.                        *)
(*                                                                      *)
(************************************************************************)


(************************************************************************)
(*                                                                      *)
(*               LAPLACIAN IN 2D, O(h^2) STENCILS                       *)
(*                                                                      *)
(************************************************************************)



(***************** General form of the stencil **************************)

Laplacian2D2hGeneral[f_,x_,y_,h_]:=(
c1(	f[x+h,y+h]+f[x+h,y-h]+f[x-h,y+h]+f[x-h,y-h])+
c2(	f[x+h,y]+f[x-h,y]+f[x,y+h]+f[x,y-h])+
c3	f[x,y]
)/h^2

cLaplacian2D2hAnisotropic={ c2->1-2 c1, c3->-4+4 c1 }

cLaplacian2D2hIsotropic={ c1->1/6, c2->2/3, c3->-10/3 }

eLaplacian2D2hIsotropic=1/12



(***************** Five-point O(h^2) anisotropic stencil ****************)

Laplacian2D2hAniso5p[f_,x_,y_,h_]:=(
f[x-h,y]+f[x+h,y]+f[x,y+h]+f[x,y-h]-4f[x,y]
)/h^2

cLaplacian2D2hAniso5p={	c1->0, c2->1, c3->-4 }




(***************** Nine-point O(h^2) isotropic stencil ******************)

Laplacian2D2hIso9p[f_,x_,y_,h_]:=(
f[x-h,y+h]+f[x+h,y-h]+f[x+h,y+h]+f[x-h,y-h]
+4f[x+h,y]+4f[x,y+h]+4f[x,y-h]+4f[x-h,y]-20f[x,y]
)/(6h^2)

cLaplacian2D2hIso9p={ c1->1/6, c2->2/3, c3->-10/3 }

eLaplacian2D2hIso9p=1/12



(************************************************************************)
(*                                                                      *)
(*               LAPLACIAN IN 2D, O(h^4) STENCILS                       *)
(*                                                                      *)
(************************************************************************)


(***************** General form of the stencil **************************)

Laplacian2D4hGeneral[f_,x_,y_,h_]:=(
c1(	f[x+2h,y+2h]+f[x-2h,y+2h]+f[x+2h,y-2h]+f[x-2h,y-2h])+
c2(	f[x+2h,y+ h]+f[x-2h,y+h]+f[x+2h,y-h]+f[x-2h,y-h]+
        f[x+ h,y+2h]+f[x-h,y+2h]+f[x+h,y-2h]+f[x-h,y-2h])+
c3(	f[x+2h,y]+f[x-2h,y]+f[x,y+2h]+f[x,y-2h])+
c4(	f[x+ h,y+h]+f[x-h,y+h]+f[x+h,y-h]+f[x-h,y-h])+
c5(	f[x+ h,y]+f[x-h,y]+f[x,y+h]+f[x,y-h])+
c6 	f[x   ,y]
)/(h^2)

cLaplacian2D4hAnisotropic={ c3->-1/12-2c1-2c2, c4->-16c1-8c2, 
c5->4/3+32c1+14c2, c6->-5-60c1-24c2 }

cLaplacian2D4hIsotropic={ c2->-1/30-4c1, c3->-1/60+6c1, c4->4/15+16c1,
c5->13/15-24c1, c6->-21/5+36c1 }

eLaplacian2D4hIsotropic= -1/90



(***************** Nine-point O(h^4) anisotropic stencil ****************)

Laplacian2D4hAniso9p[f_,x_,y_,h_]:=(
-f[x-2h,y]-f[x+2h,y]-f[x,y+2h]-f[x,y-2h]
+16f[x-h,y]+16f[x+h,y]+16f[x,y+h]+16f[x,y-h]-60f[x,y]
)/(12h^2)

cLaplacian2D4hAniso9p={ c1->0, c2->0, c4->0, c3->-1/12, c5->4/3, c6->-5 }



(***************** 17-point O(h^4) isotropic stencil ********************)

Laplacian2D4hIso17p[f_,x_,y_,h_]:=(
-1(	f[x+2h,y+2h]+f[x-2h,y+2h]+f[x+2h,y-2h]+f[x-2h,y-2h])
-8(	f[x+2h,y]+f[x-2h,y]+f[x,y+2h]+f[x,y-2h])+
16(	f[x+ h,y+h]+f[x-h,y+h]+f[x+h,y-h]+f[x-h,y-h])+
128(	f[x+ h,y]+f[x-h,y]+f[x,y+h]+f[x,y-h])
-540 	f[x   ,y]
)/(120h^2)

cLaplacian2D4hIso17p={
c1->-1/120, c2->0, c3->-1/15, c4->2/15, c5->16/15, c6->-9/2 }

eLaplacian2D4hIso17p= -1/90



(***************** 21-point O(h^4) isotropic stencil ********************)

Laplacian2D4hIso21p[f_,x_,y_,h_]:=(
-2(	f[x+2h,y+ h]+f[x-2h,y+h]+f[x+2h,y-h]+f[x-2h,y-h]+
        f[x+ h,y+2h]+f[x-h,y+2h]+f[x+h,y-2h]+f[x-h,y-2h])
-1(	f[x+2h,y]+f[x-2h,y]+f[x,y+2h]+f[x,y-2h])+
16(	f[x+ h,y+h]+f[x-h,y+h]+f[x+h,y-h]+f[x-h,y-h])+
52(	f[x+ h,y]+f[x-h,y]+f[x,y+h]+f[x,y-h])
-252 	f[x   ,y]
)/(60h^2)

cLaplacian2D4hIso21p={ 
c1->0, c2->-1/30, c3->-1/60, c4->4/15, c5->13/15, c6->-21/5 }

eLaplacian2D4hIso21p= -1/90



(************************************************************************)
(*                                                                      *)
(*               LAPLACIAN IN 3D, O(h^2) STENCILS                       *)
(*                                                                      *)
(************************************************************************)


(***************** General form of the stencil **************************)

Laplacian3D2hGeneral[f_,x_,y_,z_,h_]:=(
c1(	f[x+h,y+h,z+h]+f[x-h,y+h,z+h]+f[x+h,y-h,z+h]+f[x-h,y-h,z+h]+
        f[x+h,y+h,z-h]+f[x-h,y+h,z-h]+f[x+h,y-h,z-h]+f[x-h,y-h,z-h])+
c2(	f[x+h,y+h,z]+f[x-h,y+h,z]+f[x+h,y-h,z]+f[x-h,y-h,z]+
        f[x+h,y,z-h]+f[x-h,y,z-h]+f[x+h,y,z+h]+f[x-h,y,z+h]+
        f[x,y+h,z-h]+f[x,y+h,z+h]+f[x,y-h,z-h]+f[x,y-h,z+h])+
c3(	f[x+h,y,z]+f[x-h,y,z]+f[x,y+h,z]+f[x,y-h,z]+f[x,y,z+h]+f[x,y,z-h])+
c4	f[x,y,z]
)/(h^2)

cLaplacian3D2hAnisotropic={ c3->1-4c1-4c2, c4->-6+16c1+12c2 }

cLaplacian3D2hIsotropic={ c2->1/6-2c1, c3->1/3+4c1, c4->-4-8c1 }

eLaplacian3D2hIsotropic=1/12



(***************** Seven-point O(h^2) anisotropic stencil ***************)

Laplacian3D2hAniso7p[f_,x_,y_,z_,h_]:=(
f[x-h,y,z]+f[x+h,y,z]+f[x,y+h,z]+f[x,y-h,z]+
f[x,y,z+h]+f[x,y,z-h]-6 f[x,y,z]
)/h^2

cLaplacian3D2hAniso7p={ c1->0, c2->0, c3->1, c4->-6 }



(***************** 15-point O(h^2) isotropic stencil ********************)

Laplacian3D2hIso15p[f_,x_,y_,z_,h_]:=(
f[x-h,y-h,z-h] + f[x+h,y-h,z-h] + f[x-h,y+h,z-h] +
f[x+h,y+h,z-h] + f[x-h,y-h,z+h] + f[x+h,y-h,z+h] + f[x-h,y+h,z+h] +
f[x+h,y+h,z+h] + 
8 f[x-h,y,z]+8 f[x+h,y,z]+8 f[x,y+h,z]+8 f[x,y-h,z]+
8 f[x,y,z+h]+8 f[x,y,z-h] - 56 f[x,y,z]
)/(12 h^2)

cLaplacian3D2hIso15p={ c1->1/12, c2->0, c3->2/3, c4->-14/3 }

eLaplacian3D2hIso15p=1/12



(***************** 19-point O(h^2) isotropic stencil ********************)

Laplacian3D2hIso19p[f_,x_,y_,z_,h_]:=(
f[x-h,y-h,z]+f[x+h,y-h,z]+f[x-h,y+h,z]+f[x+h,y+h,z]+
f[x-h,y,z-h]+f[x+h,y,z-h]+f[x-h,y,z+h]+f[x+h,y,z+h]+
f[x,y-h,z-h]+f[x,y+h,z-h]+f[x,y-h,z+h]+f[x,y+h,z+h]+
2 f[x-h,y,z]+2 f[x+h,y,z]+2 f[x,y+h,z]+2 f[x,y-h,z]+
2 f[x,y,z+h]+2 f[x,y,z-h] - 24 f[x,y,z]
)/(6 h^2)

cLaplacian3D2hIso19p={ c1->0, c2->1/6, c3->1/3, c4->-4 }

eLaplacian3D2hIso19p=1/12



(***************** 21-point O(h^2) isotropic stencil ********************)

Laplacian3D2hIso21p[f_,x_,y_,z_,h_]:=(
-(f[x-h,y-h,z-h] + f[x+h,y-h,z-h] + f[x-h,y+h,z-h] +
f[x+h,y+h,z-h] + f[x-h,y-h,z+h] + f[x+h,y-h,z+h] + f[x-h,y+h,z+h] +
f[x+h,y+h,z+h]) +
4(f[x-h,y-h,z]+f[x+h,y-h,z]+f[x-h,y+h,z]+f[x+h,y+h,z]+
f[x-h,y,z-h]+f[x+h,y,z-h]+f[x-h,y,z+h]+f[x+h,y,z+h]+
f[x,y-h,z-h]+f[x,y+h,z-h]+f[x,y-h,z+h]+f[x,y+h,z+h])+
-40 f[x,y,z]
)/(12 h^2)

cLaplacian3D2hIso21p={ c1->-1/12, c2->1/3, c3->0, c4->-10/3 }

eLaplacian3D2hIso21p=1/12



(***************** 27-point O(h^2) isotropic stencil ********************)

Laplacian3D2hIso27p[f_,x_,y_,z_,h_]:=(
f[x-h,y-h,z-h] + f[x+h,y-h,z-h] + f[x-h,y+h,z-h] +
f[x+h,y+h,z-h] + f[x-h,y-h,z+h] + f[x+h,y-h,z+h] + f[x-h,y+h,z+h] +
f[x+h,y+h,z+h] +
3(f[x-h,y-h,z]+f[x+h,y-h,z]+f[x-h,y+h,z]+f[x+h,y+h,z]+
f[x-h,y,z-h]+f[x+h,y,z-h]+f[x-h,y,z+h]+f[x+h,y,z+h]+
f[x,y-h,z-h]+f[x,y+h,z-h]+f[x,y-h,z+h]+f[x,y+h,z+h])+
14( f[x-h,y,z]+ f[x+h,y,z]+ f[x,y+h,z]+ f[x,y-h,z]+
f[x,y,z+h]+ f[x,y,z-h])-128 f[x,y,z]
)/(30 h^2)

cLaplacian3D2hIso27p={ c1->1/30, c2->1/10, c3->7/15, c4->-64/15 }

eLaplacian3D2hIso27p=1/12




(************************************************************************)
(*                                                                      *)
(*               BILAPLACIAN IN 2D, O(h^2) STENCILS                     *)
(*                                                                      *)
(************************************************************************)



(***************** General form of the stencil **************************)

Bilaplacian2D2hGeneral[f_,x_,y_,h_]:=(
c1(	f[x+2h,y+2h]+f[x-2h,y+2h]+f[x+2h,y-2h]+f[x-2h,y-2h])+
c2(	f[x+2h,y+h]+f[x-2h,y+h]+f[x+2h,y-h]+f[x-2h,y-h]+
        f[x+h,y+2h]+f[x-h,y+2h]+f[x+h,y-2h]+f[x-h,y-2h])+
c3(	f[x+2h,y]+f[x-2h,y]+f[x,y+2h]+f[x,y-2h])+
c4(	f[x+h,y+h]+f[x-h,y+h]+f[x+h,y-h]+f[x-h,y-h])+
c5(	f[x+h,y]+f[x-h,y]+f[x,y+h]+f[x,y-h])+
c6	f[x,y]
)/(h^4)

cBilaplacian2D2hAnisotropic={ c3->1-2c1-2c2, c4->2-16c1-8c2, c5->32c1+14c2-8, 
c6->20-60c1-24c2 }

cBilaplacian2D2hIsotropic={ c2->1/3-4c1, c3->1/3+6c1, c4->-2/3+16c1, 
c5->-10/3-24c1, c6->12+36c1 }

eBilaplacian2D2hIsotropic=1/6



(***************** 13-point O(h^2) anisotropic stencil ******************)

Bilaplacian2D2hAniso13p[f_,x_,y_,h_]:=(
f[x,y-2h]+f[x,y+2h]+f[x-2h,y]+f[x+2h,y]+
2(f[x-h,y-h]+f[x+h,y-h]+f[x-h,y+h]+f[x+h,y+h])
-8(f[x-h,y]+f[x+h,y]+f[x,y-h]+f[x,y+h])
+20 f[x,y]
)/h^4

cBilaplacian2D2hAniso13p={ c1->0, c2->0, c3->1, c4->2, c5->-8, c6->20 }



(***************** 17-point O(h^2) isotropic stencil ********************)

Bilaplacian2D2hIso17p[f_,x_,y_,h_]:=(
(f[x+2h,y+2h]+f[x-2h,y+2h]+f[x+2h,y-2h]+f[x-2h,y-2h])
+10(f[x+2h,y]+f[x-2h,y]+f[x,y+2h]+f[x,y-2h])
+8(f[x+h,y+h]+f[x-h,y+h]+f[x+h,y-h]+f[x-h,y-h])
-64(f[x+h,y]+f[x-h,y]+f[x,y+h]+f[x,y-h])
+180 f[x,y]
)/(12 h^4)

cBilaplacian2D2hIso17p={ c1->1/12, c2->0, c3->5/6, c4->2/3, c5->-16/3, c6->15 }

eBilaplacian2D2hIso17p=1/6



(***************** 21-point O(h^2) isotropic stencil ********************)

Bilaplacian2D2hIso21p[f_,x_,y_,h_]:=(
(f[x+2h,y+h]+f[x-2h,y+h]+f[x+2h,y-h]+f[x-2h,y-h]+
        f[x+h,y+2h]+f[x-h,y+2h]+f[x+h,y-2h]+f[x-h,y-2h])
+1(f[x+2h,y]+f[x-2h,y]+f[x,y+2h]+f[x,y-2h])
-2(f[x+h,y+h]+f[x-h,y+h]+f[x+h,y-h]+f[x-h,y-h])
-10(f[x+h,y]+f[x-h,y]+f[x,y+h]+f[x,y-h])
+36 f[x,y]
)/(3 h^4)

cBilaplacian2D2hIso21p={ c1->0, c2->1/3, c3->1/3, c4->-2/3, c5->-10/3, c6->12 }

eBilaplacian2D2hIso21p=1/6



(************************************************************************)
(*                                                                      *)
(*               BILAPLACIAN IN 2D, O(h^4) STENCILS                     *)
(*                                                                      *)
(************************************************************************)



(***************** General form of the stencil **************************)

Bilaplacian2D4hGeneral[f_,x_,y_,h_]:=(
c1(	f[x+3h,y+3h]+f[x-3h,y+3h]+f[x+3h,y-3h]+f[x-3h,y-3h])+
c2(	f[x+3h,y+2h]+f[x-3h,y+2h]+f[x+3h,y-2h]+f[x-3h,y-2h]+
        f[x+2h,y+3h]+f[x-2h,y+3h]+f[x+2h,y-3h]+f[x-2h,y-3h])+
c3(	f[x+3h,y+h]+f[x-3h,y+h]+f[x+3h,y-h]+f[x-3h,y-h]+
        f[x+h,y+3h]+f[x-h,y+3h]+f[x+h,y-3h]+f[x-h,y-3h])+
c4(	f[x+3h,y]+f[x-3h,y]+f[x,y+3h]+f[x,y-3h])+
c5(	f[x+2h,y+2h]+f[x-2h,y+2h]+f[x+2h,y-2h]+f[x-2h,y-2h])+
c6(	f[x+2h,y+h]+f[x-2h,y+h]+f[x+2h,y-h]+f[x-2h,y-h]+
        f[x+h,y+2h]+f[x-h,y+2h]+f[x+h,y-2h]+f[x-h,y-2h])+
c7(	f[x+2h,y]+f[x-2h,y]+f[x,y+2h]+f[x,y-2h])+
c8(	f[x+h,y+h]+f[x-h,y+h]+f[x+h,y-h]+f[x-h,y-h])+
c9(	f[x+h,y]+f[x-h,y]+f[x,y+h]+f[x,y-h])+
c10	f[x,y]
)/(h^4)

cBilaplacian2D4hAnisotropic={ c4->-1/6-2c1-2c2-2c3, c6->-1/6-54c1-33c2-6c3-4c5,
c7->7/3+108c1+64c2+12c3+6c5, c8->10/3+351c1+192c2+30c3+16c5,
c9->-77/6-594c1-318c2-50c3-24c5, c10->92/3+976c1+512c2+80c3+36c5 }

cBilaplacian2D4hIsotropic={ c3->-17/180-9c1-4c2, c4->1/45+16c1+6c2,
c5->-29/180-36c1-12c2, c6->47/45+144c1+39c2, c7->7/30-216c1-56c2,
c8->-187/90-495c1-120c2, c9->-191/45+720c1+170c2, c10->779/45-1040c1-240c2 }

eBilaplacian2D4hIsotropic= -7/240




(***************** 21-point O(h^4) anisotropic stencil ******************)

Bilaplacian2D4hAniso21p[f_,x_,y_,h_]:=(
-4(	f[x+3h,y]+f[x-3h,y]+f[x,y+3h]+f[x,y-3h])
-1(	f[x+2h,y+2h]+f[x-2h,y+2h]+f[x+2h,y-2h]+f[x-2h,y-2h])
+50(	f[x+2h,y]+f[x-2h,y]+f[x,y+2h]+f[x,y-2h])
+64(	f[x+h,y+h]+f[x-h,y+h]+f[x+h,y-h]+f[x-h,y-h])
-284(	f[x+h,y]+f[x-h,y]+f[x,y+h]+f[x,y-h])
+700	f[x,y]
)/(24 h^4)

cBilaplacian2D4hAniso21p={ c1->0, c2->0, c3->0, c4->-1/6, c5->-1/24, 
c6->0, c7->25/12, c8->8/3, c9->-71/6, c10->175/6 }



(***************** 25-point O(h^4) anisotropic stencil ******************)

Bilaplacian2D4hAniso25p[f_,x_,y_,h_]:=(
-(	f[x+3h,y]+f[x-3h,y]+f[x,y+3h]+f[x,y-3h])
-(	f[x+2h,y+h]+f[x-2h,y+h]+f[x+2h,y-h]+f[x-2h,y-h]+
        f[x+h,y+2h]+f[x-h,y+2h]+f[x+h,y-2h]+f[x-h,y-2h])+
14(	f[x+2h,y]+f[x-2h,y]+f[x,y+2h]+f[x,y-2h])+
20(	f[x+h,y+h]+f[x-h,y+h]+f[x+h,y-h]+f[x-h,y-h])
-77(	f[x+h,y]+f[x-h,y]+f[x,y+h]+f[x,y-h])+
184	f[x,y]
)/(6 h^4)

cBilaplacian2D4hAniso25p={ c1->0, c2->0, c3->0, c4->-1/6, c5->0, c6->-1/6,
c7->7/3, c8->10/3, c9->-77/6, c10->92/3 }



(***************** 33-point O(h^4) isotropic stencil ********************)

Bilaplacian2D4hIso33p[f_,x_,y_,h_]:=(
-17(f[x+3h,y+3h]+f[x-3h,y+3h]+f[x+3h,y-3h]+f[x-3h,y-3h])
-236(f[x+3h,y]+f[x-3h,y]+f[x,y+3h]+f[x,y-3h])
+351(f[x+2h,y+2h]+f[x-2h,y+2h]+f[x+2h,y-2h]+f[x-2h,y-2h])
-756(f[x+2h,y+h]+f[x-2h,y+h]+f[x+2h,y-h]+f[x-2h,y-h]+
        f[x+h,y+2h]+f[x-h,y+2h]+f[x+h,y-2h]+f[x-h,y-2h])
+4050(f[x+2h,y]+f[x-2h,y]+f[x,y+2h]+f[x,y-2h])
+5049(f[x+h,y+h]+f[x-h,y+h]+f[x+h,y-h]+f[x-h,y-h])
-19116(f[x+h,y]+f[x-h,y]+f[x,y+h]+f[x,y-h])
+45724 f[x,y]
)/(1620 h^4)

cBilaplacian2D4hIso33p={ c1->-17/1620, c2->0, c3->0, c4->-59/405, c5->13/60,
c6->-7/15, c7->5/2, c8->187/60, c9->-59/5, c10->11431/405 }

eBilaplacian2D4hIso33p= -7/240



(***************** 37-point O(h^4) isotropic stencil ********************)

Bilaplacian2D4hIso37p[f_,x_,y_,h_]:=(
-17(f[x+3h,y+h]+f[x-3h,y+h]+f[x+3h,y-h]+f[x-3h,y-h]+
        f[x+h,y+3h]+f[x-h,y+3h]+f[x+h,y-3h]+f[x-h,y-3h])
+4(f[x+3h,y]+f[x-3h,y]+f[x,y+3h]+f[x,y-3h])
-29(f[x+2h,y+2h]+f[x-2h,y+2h]+f[x+2h,y-2h]+f[x-2h,y-2h])
+188(f[x+2h,y+h]+f[x-2h,y+h]+f[x+2h,y-h]+f[x-2h,y-h]+
        f[x+h,y+2h]+f[x-h,y+2h]+f[x+h,y-2h]+f[x-h,y-2h])
+42(f[x+2h,y]+f[x-2h,y]+f[x,y+2h]+f[x,y-2h])
-374(f[x+h,y+h]+f[x-h,y+h]+f[x+h,y-h]+f[x-h,y-h])
-764(f[x+h,y]+f[x-h,y]+f[x,y+h]+f[x,y-h])
+3116 f[x,y]
)/(180 h^4)

cBilaplacian2D4hIso37p={ c1->0, c2->0, c3->-17/180, c4->1/45, c5->-29/180,
c6->47/45, c7->7/30, c8->-187/90, c9->-191/45, c10->779/45 }

eBilaplacian2D4hIso37p= -7/240



(************************************************************************)
(*                                                                      *)
(*               BILAPLACIAN IN 3D, O(h^2) STENCILS                     *)
(*                                                                      *)
(************************************************************************)



(***************** General form of the stencil **************************)

Bilaplacian3D2hGeneral[f_,x_,y_,z_,h_]:=(
c1(	f[x+2h,y+2h,z+2h]+f[x-2h,y+2h,z+2h]+f[x+2h,y-2h,z+2h]+f[x-2h,y-2h,z+2h]+
	f[x+2h,y+2h,z-2h]+f[x-2h,y+2h,z-2h]+f[x+2h,y-2h,z-2h]+f[x-2h,y-2h,z-2h])
+c2(	f[x+2h,y+2h,z+h]+f[x-2h,y+2h,z+h]+f[x+2h,y-2h,z+h]+f[x-2h,y-2h,z+h]+
        f[x+2h,y+2h,z-h]+f[x-2h,y+2h,z-h]+f[x+2h,y-2h,z-h]+f[x-2h,y-2h,z-h]+
        f[x+2h,y+h,z+2h]+f[x-2h,y+h,z+2h]+f[x+2h,y-h,z+2h]+f[x-2h,y-h,z+2h]+
        f[x+2h,y+h,z-2h]+f[x-2h,y+h,z-2h]+f[x+2h,y-h,z-2h]+f[x-2h,y-h,z-2h]+
        f[x+h,y+2h,z+2h]+f[x-h,y+2h,z+2h]+f[x+h,y-2h,z+2h]+f[x-h,y-2h,z+2h]+
        f[x+h,y+2h,z-2h]+f[x-h,y+2h,z-2h]+f[x+h,y-2h,z-2h]+f[x-h,y-2h,z-2h])+
+c3(	f[x+2h,y+2h,z]+f[x-2h,y+2h,z]+f[x+2h,y-2h,z]+f[x-2h,y-2h,z]+
        f[x+2h,y,z-2h]+f[x-2h,y,z-2h]+f[x+2h,y,z+2h]+f[x-2h,y,z+2h]+
        f[x,y+2h,z-2h]+f[x,y+2h,z+2h]+f[x,y-2h,z-2h]+f[x,y-2h,z+2h])
+c4(	f[x+h,y+h,z+2h]+f[x-h,y+h,z+2h]+f[x+h,y-h,z+2h]+f[x-h,y-h,z+2h]+
        f[x+h,y+h,z-2h]+f[x-h,y+h,z-2h]+f[x+h,y-h,z-2h]+f[x-h,y-h,z-2h]+
        f[x+h,y+2h,z+h]+f[x-h,y+2h,z+h]+f[x+h,y-2h,z+h]+f[x-h,y-2h,z+h]+
        f[x+h,y+2h,z-h]+f[x-h,y+2h,z-h]+f[x+h,y-2h,z-h]+f[x-h,y-2h,z-h]+
        f[x+2h,y+h,z+h]+f[x-2h,y+h,z+h]+f[x+2h,y-h,z+h]+f[x-2h,y-h,z+h]+
        f[x+2h,y+h,z-h]+f[x-2h,y+h,z-h]+f[x+2h,y-h,z-h]+f[x-2h,y-h,z-h])
+c5(	f[x,y+h,z+2h]+f[x,y-h,z+2h]+f[x,y+h,z-2h]+f[x,y-h,z-2h]+
        f[x,y+2h,z+h]+f[x,y-2h,z+h]+f[x,y+2h,z-h]+f[x,y-2h,z-h]+
        f[x+2h,y,z+h]+f[x-2h,y,z+h]+f[x+2h,y,z-h]+f[x-2h,y,z-h]+
        f[x+h,y,z+2h]+f[x-h,y,z+2h]+f[x+h,y,z-2h]+f[x-h,y,z-2h]+
        f[x+h,y+2h,z]+f[x-h,y+2h,z]+f[x+h,y-2h,z]+f[x-h,y-2h,z]+
        f[x+2h,y+h,z]+f[x-2h,y+h,z]+f[x+2h,y-h,z]+f[x-2h,y-h,z])
+c6(	f[x+h,y+h,z+h]+f[x-h,y+h,z+h]+f[x+h,y-h,z+h]+f[x-h,y-h,z+h]+
        f[x+h,y+h,z-h]+f[x-h,y+h,z-h]+f[x+h,y-h,z-h]+f[x-h,y-h,z-h])
+c7(	f[x+h,y+h,z]+f[x-h,y+h,z]+f[x+h,y-h,z]+f[x-h,y-h,z]+
        f[x+h,y,z-h]+f[x-h,y,z-h]+f[x+h,y,z+h]+f[x-h,y,z+h]+
        f[x,y+h,z-h]+f[x,y+h,z+h]+f[x,y-h,z-h]+f[x,y-h,z+h])
+c8(	f[x+2h,y,z]+f[x-2h,y,z]+f[x,y+2h,z]+f[x,y-2h,z]+f[x,y,z+2h]+f[x,y,z-2h])
+c9(	f[x+h,y,z]+f[x-h,y,z]+f[x,y+h,z]+f[x,y-h,z]+f[x,y,z+h]+f[x,y,z-h])
+c10	f[x,y,z]
)/(h^4)

cBilaplacian3D2hAnisotropic={ c7->2-32c1-48c2-16c3-18c4-8c5-2c6,
c8->1-4c1-8c2-4c3-4c4-4c5, c9->128c1+188c2+64c3+64c4+28c5+4c6-12,
c10->42-368c1-528c2-180c3-168c4-72c5-8c6 }

cBilaplacian3D2hIsotropic={ c1->1/24(1-30c2-12c3-6c4-3c5),
c6->32c2+32c3+4c4+8c5-5/3, c7->4-72c2-64c3-18c4-20c5, c8->5/6-3c2-2c3-3c4-7/2c5,
c9->-40/3+156c2+128c3+48c4+44c5, c10->40-324c2-252c3-108c4-90c5 }

eBilaplacian3D2hIsotropic=1/6



(***************** 21-point O(h^2) anisotropic stencil ******************)

Bilaplacian3D2hAniso21p[f_,x_,y_,z_,h_]:=(
(	f[x+2h,y+2h,z]+f[x-2h,y+2h,z]+f[x+2h,y-2h,z]+f[x-2h,y-2h,z]+
        f[x+2h,y,z-2h]+f[x-2h,y,z-2h]+f[x+2h,y,z+2h]+f[x-2h,y,z+2h]+
        f[x,y+2h,z-2h]+f[x,y+2h,z+2h]+f[x,y-2h,z-2h]+f[x,y-2h,z+2h])
-4(	f[x+h,y+h,z+h]+f[x-h,y+h,z+h]+f[x+h,y-h,z+h]+f[x-h,y-h,z+h]+
        f[x+h,y+h,z-h]+f[x-h,y+h,z-h]+f[x+h,y-h,z-h]+f[x-h,y-h,z-h])
+20	f[x,y,z]
)/(4 h^4)

cBilaplacian3D2hAniso21p={ c1->0, c2->0, c3->1/4, c4->0, c5->0, c6->-1, 
c7->0, c8->0, c9->0, c10->5 }



(***************** 25-point O(h^2) anisotropic stencil ******************)

Bilaplacian3D2hAniso25p[f_,x_,y_,z_,h_]:=(
2(	f[x+h,y+h,z]+f[x-h,y+h,z]+f[x+h,y-h,z]+f[x-h,y-h,z]+
	f[x+h,y,z-h]+f[x-h,y,z-h]+f[x+h,y,z+h]+f[x-h,y,z+h]+
	f[x,y+h,z-h]+f[x,y+h,z+h]+f[x,y-h,z-h]+f[x,y-h,z+h])
+(	f[x+2h,y,z]+f[x-2h,y,z]+f[x,y+2h,z]+f[x,y-2h,z]+f[x,y,z+2h]+f[x,y,z-2h])
-12(	f[x+h,y,z]+f[x-h,y,z]+f[x,y+h,z]+f[x,y-h,z]+f[x,y,z+h]+f[x,y,z-h])
+42	f[x,y,z]
)/(h^4)

cBilaplacian3D2hAniso25p={ c1->0, c2->0, c3->0, c4->0, c5->0, c6->0,
c7->2, c8->1, c9->-12, c10->42 }




(***************** 41-point O(h^2) isotropic stencil ********************)

Bilaplacian3D2hIso41p[f_,x_,y_,z_,h_]:=(
(	f[x+2h,y+2h,z+2h]+f[x-2h,y+2h,z+2h]+f[x+2h,y-2h,z+2h]+f[x-2h,y-2h,z+2h]+
        f[x+2h,y+2h,z-2h]+f[x-2h,y+2h,z-2h]+f[x+2h,y-2h,z-2h]+f[x-2h,y-2h,z-2h])
-40(	f[x+h,y+h,z+h]+f[x-h,y+h,z+h]+f[x+h,y-h,z+h]+f[x-h,y-h,z+h]+
        f[x+h,y+h,z-h]+f[x-h,y+h,z-h]+f[x+h,y-h,z-h]+f[x-h,y-h,z-h])
+96(	f[x+h,y+h,z]+f[x-h,y+h,z]+f[x+h,y-h,z]+f[x-h,y-h,z]+
        f[x+h,y,z-h]+f[x-h,y,z-h]+f[x+h,y,z+h]+f[x-h,y,z+h]+
        f[x,y+h,z-h]+f[x,y+h,z+h]+f[x,y-h,z-h]+f[x,y-h,z+h])
+20(	f[x+2h,y,z]+f[x-2h,y,z]+f[x,y+2h,z]+f[x,y-2h,z]+f[x,y,z+2h]+f[x,y,z-2h])
-320(	f[x+h,y,z]+f[x-h,y,z]+f[x,y+h,z]+f[x,y-h,z]+f[x,y,z+h]+f[x,y,z-h])
+960	f[x,y,z]
)/(24 h^4)

cBilaplacian3D2hIso41p={ c1->1/24, c2->0, c3->0, c4->0, c5->0, c6->-5/3,
c7->4, c8->5/6, c9->-40/3, c10->40 }

eBilaplacian3D2hIso41p=1/6



(***************** 52-point O(h^2) isotropic stencil ********************)

Bilaplacian3D2hIso52p[f_,x_,y_,z_,h_]:=(
-1(	f[x+2h,y+2h,z+2h]+f[x-2h,y+2h,z+2h]+f[x+2h,y-2h,z+2h]+f[x-2h,y-2h,z+2h]+
	f[x+2h,y+2h,z-2h]+f[x-2h,y+2h,z-2h]+f[x+2h,y-2h,z-2h]+f[x-2h,y-2h,z-2h])
+10(	f[x+h,y+h,z+2h]+f[x-h,y+h,z+2h]+f[x+h,y-h,z+2h]+f[x-h,y-h,z+2h]+
        f[x+h,y+h,z-2h]+f[x-h,y+h,z-2h]+f[x+h,y-h,z-2h]+f[x-h,y-h,z-2h]+
        f[x+h,y+2h,z+h]+f[x-h,y+2h,z+h]+f[x+h,y-2h,z+h]+f[x-h,y-2h,z+h]+
        f[x+h,y+2h,z-h]+f[x-h,y+2h,z-h]+f[x+h,y-2h,z-h]+f[x-h,y-2h,z-h]+
        f[x+2h,y+h,z+h]+f[x-2h,y+h,z+h]+f[x+2h,y-h,z+h]+f[x-2h,y-h,z+h]+
        f[x+2h,y+h,z-h]+f[x-2h,y+h,z-h]+f[x+2h,y-h,z-h]+f[x-2h,y-h,z-h])
-20(	f[x+h,y+h,z+h]+f[x-h,y+h,z+h]+f[x+h,y-h,z+h]+f[x-h,y-h,z+h]+
        f[x+h,y+h,z-h]+f[x-h,y+h,z-h]+f[x+h,y-h,z-h]+f[x-h,y-h,z-h])
-36(	f[x+h,y+h,z]+f[x-h,y+h,z]+f[x+h,y-h,z]+f[x-h,y-h,z]+
        f[x+h,y,z-h]+f[x-h,y,z-h]+f[x+h,y,z+h]+f[x-h,y,z+h]+
        f[x,y+h,z-h]+f[x,y+h,z+h]+f[x,y-h,z-h]+f[x,y-h,z+h])
+360	f[x,y,z]
)/(36 h^4)

cBilaplacian3D2hIso52p={ c1->-1/36, c2->0, c3->0, c4->5/18, c5->0, c6->-5/9, 
c7->-1, c8->0, c9->0, c10->10 }

eBilaplacian3D2hIso52p=1/6



(***************** 57-point O(h^2) isotropic stencil ********************)

Bilaplacian3D2hIso57p[f_,x_,y_,z_,h_]:=(
(	f[x,y+h,z+2h]+f[x,y-h,z+2h]+f[x,y+h,z-2h]+f[x,y-h,z-2h]+
        f[x,y+2h,z+h]+f[x,y-2h,z+h]+f[x,y+2h,z-h]+f[x,y-2h,z-h]+
        f[x+2h,y,z+h]+f[x-2h,y,z+h]+f[x+2h,y,z-h]+f[x-2h,y,z-h]+
        f[x+h,y,z+2h]+f[x-h,y,z+2h]+f[x+h,y,z-2h]+f[x-h,y,z-2h]+
        f[x+h,y+2h,z]+f[x-h,y+2h,z]+f[x+h,y-2h,z]+f[x-h,y-2h,z]+
        f[x+2h,y+h,z]+f[x-2h,y+h,z]+f[x+2h,y-h,z]+f[x-2h,y-h,z])
+3(	f[x+h,y+h,z+h]+f[x-h,y+h,z+h]+f[x+h,y-h,z+h]+f[x-h,y-h,z+h]+
        f[x+h,y+h,z-h]+f[x-h,y+h,z-h]+f[x+h,y-h,z-h]+f[x-h,y-h,z-h])
-8(	f[x+h,y+h,z]+f[x-h,y+h,z]+f[x+h,y-h,z]+f[x-h,y-h,z]+
        f[x+h,y,z-h]+f[x-h,y,z-h]+f[x+h,y,z+h]+f[x-h,y,z+h]+
        f[x,y+h,z-h]+f[x,y+h,z+h]+f[x,y-h,z-h]+f[x,y-h,z+h])
-1(	f[x+2h,y,z]+f[x-2h,y,z]+f[x,y+2h,z]+f[x,y-2h,z]+f[x,y,z+2h]+f[x,y,z-2h])
+4(	f[x+h,y,z]+f[x-h,y,z]+f[x,y+h,z]+f[x,y-h,z]+f[x,y,z+h]+f[x,y,z-h])
+30	f[x,y,z]
)/(3 h^4)

cBilaplacian3D2hIso57p={ c1->0, c2->0, c3->0, c4->0, c5->1/3, c6->1, c7->-8/3,
c8->-1/3, c9->4/3, c10->10 }

eBilaplacian3D2hIso57p=1/6





(************************************************************************)
(*                                                                      *)
(*        GRADIENT OF A LAPLACIAN IN 2D, O(h^2) STENCILS                *)
(*                                                                      *)
(************************************************************************)



(***************** General form of the stencil **************************)

GradLap2D2hGeneral[f_,x_,y_,h_]:=(
+c1(	f[x+2h,y+2h]-f[x-2h,y+2h]+f[x+2h,y-2h]-f[x-2h,y-2h])
+c2(	f[x+2h,y+h]-f[x-2h,y+h]+f[x+2h,y-h]-f[x-2h,y-h])
+c3(	f[x+2h,y]-f[x-2h,y])
+c4(	f[x+h,y+2h]-f[x-h,y+2h]+f[x+h,y-2h]-f[x-h,y-2h])
+c5(	f[x+h,y+h]-f[x-h,y+h]+f[x+h,y-h]-f[x-h,y-h])
+c6(	f[x+h,y]-f[x-h,y])
)/h^3

cGradLap2D2hAnisotropic={ c3->1/2-2c1-2c2, c5->1/2-8c1-2c2-4c4,
c6->-2+16c1+4c2+6c4 }

cGradLap2D2hIsotropic={ c2->1/6-4c1, c3->1/6+6c1, c4->1/12-2c1,
c5->-1/6+8c1, c6->-5/6-12c1 }

eGradLap2D2hIsotropic=1/4



(***************** Six-point O(h^2) anisotropic stencil *****************)

GradLap2D2hAniso6p[f_,x_,y_,h_]:=(
(	f[x+2h,y+h]-f[x-2h,y+h]+f[x+2h,y-h]-f[x-2h,y-h])
-4(	f[x+h,y]-f[x-h,y])
)/(4h^3)

cGradLap2D2hAniso6p={ c1->0, c2->1/4, c3->0, c4->0, c5->0, c6->-1 }



(***************** Eight-point O(h^2) anisotropic stencil ***************)

GradLap2D2hAniso8p[f_,x_,y_,h_]:=(
(	f[x+2h,y]-f[x-2h,y])
+(	f[x+h,y+h]-f[x-h,y+h]+f[x+h,y-h]-f[x-h,y-h])
-4(	f[x+h,y]-f[x-h,y])
)/(2h^3)

cGradLap2D2hAniso8p={ c1->0, c2->0, c3->1/2, c4->0, c5->1/2, c6->-2 }



(***************** 12-point O(h^2) isotropic stencil ********************)

GradLap2D2hIso12p[f_,x_,y_,h_]:=(
(	f[x+2h,y+2h]-f[x-2h,y+2h]+f[x+2h,y-2h]-f[x-2h,y-2h])
+10(	f[x+2h,y]-f[x-2h,y])
+4(	f[x+h,y+h]-f[x-h,y+h]+f[x+h,y-h]-f[x-h,y-h])
-32(	f[x+h,y]-f[x-h,y])
)/(24 h^3)

cGradLap2D2hIso12p={ c1->1/24, c2->0, c3->5/12, c4->0, c5->1/6, c6->-4/3 }

eGradLap2D2hIso12p=1/4



(***************** 16-point O(h^2) isotropic stencil ********************)

GradLap2D2hIso16p[f_,x_,y_,h_]:=(
2(	f[x+2h,y+h]-f[x-2h,y+h]+f[x+2h,y-h]-f[x-2h,y-h])
+2(	f[x+2h,y]-f[x-2h,y])
+(	f[x+h,y+2h]-f[x-h,y+2h]+f[x+h,y-2h]-f[x-h,y-2h])
-2(	f[x+h,y+h]-f[x-h,y+h]+f[x+h,y-h]-f[x-h,y-h])
-10(	f[x+h,y]-f[x-h,y])
)/(12 h^3)

cGradLap2D2hIso16p={ c1->0, c2->1/6, c3->1/6, c4->1/12, c5->-1/6, c6->-5/6}

eGradLap2D2hIso16p=1/4




(************************************************************************)
(*                                                                      *)
(*        GRADIENT OF A LAPLACIAN IN 2D, O(h^4) STENCILS                *)
(*                                                                      *)
(************************************************************************)


(***************** General form of the stencil **************************)

GradLap2D4hGeneral[f_,x_,y_,h_]:=(
+c1(	f[x+3h,y+3h]-f[x-3h,y+3h]+f[x+3h,y-3h]-f[x-3h,y-3h])
+c2(	f[x+3h,y+2h]-f[x-3h,y+2h]+f[x+3h,y-2h]-f[x-3h,y-2h])
+c3(	f[x+3h,y +h]-f[x-3h,y +h]+f[x+3h,y -h]-f[x-3h,y -h])
+c4(	f[x+3h,y   ]-f[x-3h,y   ])
+c5(	f[x+2h,y+3h]-f[x-2h,y+3h]+f[x+2h,y-3h]-f[x-2h,y-3h])
+c6(	f[x+2h,y+2h]-f[x-2h,y+2h]+f[x+2h,y-2h]-f[x-2h,y-2h])
+c7(	f[x+2h,y +h]-f[x-2h,y +h]+f[x+2h,y -h]-f[x-2h,y -h])
+c8(	f[x+2h,y   ]-f[x-2h,y   ])
+c9(	f[x+ h,y+3h]-f[x- h,y+3h]+f[x+ h,y-3h]-f[x- h,y-3h])
+c10(	f[x+ h,y+2h]-f[x- h,y+2h]+f[x+ h,y-2h]-f[x- h,y-2h])
+c11(	f[x+ h,y +h]-f[x- h,y +h]+f[x+ h,y -h]-f[x- h,y -h])
+c12(	f[x+ h,y   ]-f[x- h,y   ])
)/h^3

cGradLap2D4hAnisotropic={ c4->-1/8-2c1-2c2-2c3, c7->-1/12-36c1-16c2-4c3-9c5-4c6,
c8->7/6+72c1+32c2+8c3+16c5+6c6, c10->-1/24-18c1-3c2-12c5-2c6-6c9,
c11->5/6+117c1+32c2+5c3+48c5+8c6+15c9,
c12->-77/24-198c1-58c2-10c3-72c5-12c6-20c9 }

cGradLap2D4hIsotropic={ c3->-17/240-9c1-4c2, c4->1/60+16c1+6c2,
c6->-29/360-24c1-4c2-6c5, c7->47/90+96c1+16c2+15c5,
c8->7/60-144c1-24c2-20c5, c9->-17/720-3c1-2c5, c10->47/180+48c1+5c2+12c5,
c11->-187/360-165c1-20c2-30c5, c12->-191/180+240c1+30c2+40c5 }

eGradLap2D4hIsotropic=-7/120




(***************** 14-point O(h^4) anisotropic stencil ******************)

GradLap2D4hAniso14p[f_,x_,y_,h_]:=(
-6(	f[x+3h,y   ]-f[x-3h,y   ])
-1(	f[x+2h,y+2h]-f[x-2h,y+2h]+f[x+2h,y-2h]-f[x-2h,y-2h])
+50(	f[x+2h,y   ]-f[x-2h,y   ])
+32(	f[x+ h,y +h]-f[x- h,y +h]+f[x+ h,y -h]-f[x- h,y -h])
-142(	f[x+ h,y   ]-f[x- h,y   ])
)/(48 h^3)

cGradLap2D4hAniso14p={ c1->0, c2->0, c3->0, c4->-1/8, c5->0, c6->-1/48, 
c7->0, c8->25/24, c9->0, c10->0, c11->2/3, c12->-71/24 }
 


(***************** 18-point O(h^4) anisotropic stencil ******************)

GradLap2D4hAniso18p[f_,x_,y_,h_]:=(
-3(	f[x+3h,y   ]-f[x-3h,y   ])
-2(	f[x+2h,y +h]-f[x-2h,y +h]+f[x+2h,y -h]-f[x-2h,y -h])
+28(	f[x+2h,y   ]-f[x-2h,y   ])
-1(	f[x+ h,y+2h]-f[x- h,y+2h]+f[x+ h,y-2h]-f[x- h,y-2h])
+20(	f[x+ h,y +h]-f[x- h,y +h]+f[x+ h,y -h]-f[x- h,y -h])
-77(	f[x+ h,y   ]-f[x- h,y   ])
)/(24 h^3)

cGradLap2D4hAniso18p={ c1->0, c2->0, c3->0, c4->-1/8, c5->0, c6->0, c7->-1/12,
c8->7/6, c9->0, c10->-1/24, c11->5/6, c12->-77/24 }



(***************** 26-point O(h^4) isotropic stencil ********************)

GradLap2D4hIso26p[f_,x_,y_,h_]:=(
-17(	f[x+3h,y+3h]-f[x-3h,y+3h]+f[x+3h,y-3h]-f[x-3h,y-3h])
-236(	f[x+3h,y   ]-f[x-3h,y   ])
+234(	f[x+2h,y+2h]-f[x-2h,y+2h]+f[x+2h,y-2h]-f[x-2h,y-2h])
-504(	f[x+2h,y +h]-f[x-2h,y +h]+f[x+2h,y -h]-f[x-2h,y -h])
+2700(	f[x+2h,y   ]-f[x-2h,y   ])
-252(	f[x+ h,y+2h]-f[x- h,y+2h]+f[x+ h,y-2h]-f[x- h,y-2h])
+1683(	f[x+ h,y +h]-f[x- h,y +h]+f[x+ h,y -h]-f[x- h,y -h])
-6372(	f[x+ h,y   ]-f[x- h,y   ])
)/(2160 h^3)

cGradLap2D4hIso26p={ c1->-17/2160, c2->0, c3->0, c4->-59/540, c5->0, c6->13/120,
c7->-7/30, c8->5/4, c9->0, c10->-7/60, c11->187/240, c12->-59/20 }

eGradLap2D4hIso26p=-7/120



(***************** 30-point O(h^4) isotropic stencil ********************)

GradLap2D4hIso30p[f_,x_,y_,h_]:=(
-51(	f[x+3h,y +h]-f[x-3h,y +h]+f[x+3h,y -h]-f[x-3h,y -h])
+12(	f[x+3h,y   ]-f[x-3h,y   ])
-58(	f[x+2h,y+2h]-f[x-2h,y+2h]+f[x+2h,y-2h]-f[x-2h,y-2h])
+376(	f[x+2h,y +h]-f[x-2h,y +h]+f[x+2h,y -h]-f[x-2h,y -h])
+84(	f[x+2h,y   ]-f[x-2h,y   ])
-17(	f[x+ h,y+3h]-f[x- h,y+3h]+f[x+ h,y-3h]-f[x- h,y-3h])
+188(	f[x+ h,y+2h]-f[x- h,y+2h]+f[x+ h,y-2h]-f[x- h,y-2h])
-374(	f[x+ h,y +h]-f[x- h,y +h]+f[x+ h,y -h]-f[x- h,y -h])
-764(	f[x+ h,y   ]-f[x- h,y   ])
)/(720 h^3)

cGradLap2D4hIso30p={ c1->0, c2->0, c3->-17/240, c4->1/60, c5->0, c6->-29/360,
c7->47/90, c8->7/60, c9->-17/720, c10->47/180, c11->-187/360, c12->-191/180 }

eGradLap2D4hIso30p=-7/120





(************************************************************************)
(*                                                                      *)
(*        GRADIENT OF A LAPLACIAN IN 3D, O(h^2) STENCILS                *)
(*                                                                      *)
(************************************************************************)


(***************** General form of the stencil **************************)

GradLap3D2hGeneral[f_,x_,y_,z_,h_]:=(
c1(     f[x+2h,y+2h,z+2h]-f[x-2h,y+2h,z+2h]+f[x+2h,y-2h,z+2h]-f[x-2h,y-2h,z+2h]+
        f[x+2h,y+2h,z-2h]-f[x-2h,y+2h,z-2h]+f[x+2h,y-2h,z-2h]-f[x-2h,y-2h,z-2h])
+c2(	f[x+2h,y+2h,z+h]-f[x-2h,y+2h,z+h]+f[x+2h,y-2h,z+h]-f[x-2h,y-2h,z+h]+
        f[x+2h,y+2h,z-h]-f[x-2h,y+2h,z-h]+f[x+2h,y-2h,z-h]-f[x-2h,y-2h,z-h]+
        f[x+2h,y+h,z+2h]-f[x-2h,y+h,z+2h]+f[x+2h,y-h,z+2h]-f[x-2h,y-h,z+2h]+
        f[x+2h,y+h,z-2h]-f[x-2h,y+h,z-2h]+f[x+2h,y-h,z-2h]-f[x-2h,y-h,z-2h])
+c3(    f[x+2h,y+h,z+h]-f[x-2h,y+h,z+h]+f[x+2h,y-h,z+h]-f[x-2h,y-h,z+h]+
        f[x+2h,y+h,z-h]-f[x-2h,y+h,z-h]+f[x+2h,y-h,z-h]-f[x-2h,y-h,z-h])
+c4(    f[x+2h,y+2h,z]-f[x-2h,y+2h,z]+f[x+2h,y-2h,z]-f[x-2h,y-2h,z]+
        f[x+2h,y,z-2h]-f[x-2h,y,z-2h]+f[x+2h,y,z+2h]-f[x-2h,y,z+2h])
+c5(    f[x+2h,y,z+h]-f[x-2h,y,z+h]+f[x+2h,y,z-h]-f[x-2h,y,z-h]+
        f[x+2h,y+h,z]-f[x-2h,y+h,z]+f[x+2h,y-h,z]-f[x-2h,y-h,z])
+c6(    f[x+2h,y,z]-f[x-2h,y,z])
+c7(    f[x+h,y+2h,z+2h]-f[x-h,y+2h,z+2h]+f[x+h,y-2h,z+2h]-f[x-h,y-2h,z+2h]+
        f[x+h,y+2h,z-2h]-f[x-h,y+2h,z-2h]+f[x+h,y-2h,z-2h]-f[x-h,y-2h,z-2h])
+c8(	f[x+h,y+2h,z+h]-f[x-h,y+2h,z+h]+f[x+h,y-2h,z+h]-f[x-h,y-2h,z+h]+
        f[x+h,y+2h,z-h]-f[x-h,y+2h,z-h]+f[x+h,y-2h,z-h]-f[x-h,y-2h,z-h]+
	f[x+h,y+h,z+2h]-f[x-h,y+h,z+2h]+f[x+h,y-h,z+2h]-f[x-h,y-h,z+2h]+
        f[x+h,y+h,z-2h]-f[x-h,y+h,z-2h]+f[x+h,y-h,z-2h]-f[x-h,y-h,z-2h])
+c9(    f[x+h,y,z+2h]-f[x-h,y,z+2h]+f[x+h,y,z-2h]-f[x-h,y,z-2h]+
        f[x+h,y+2h,z]-f[x-h,y+2h,z]+f[x+h,y-2h,z]-f[x-h,y-2h,z])
+c10(   f[x+h,y+h,z+h]-f[x-h,y+h,z+h]+f[x+h,y-h,z+h]-f[x-h,y-h,z+h]+
        f[x+h,y+h,z-h]-f[x-h,y+h,z-h]+f[x+h,y-h,z-h]-f[x-h,y-h,z-h])
+c11(   f[x+h,y+h,z]-f[x-h,y+h,z]+f[x+h,y-h,z]-f[x-h,y-h,z]+
        f[x+h,y,z-h]-f[x-h,y,z-h]+f[x+h,y,z+h]-f[x-h,y,z+h])
+c12(   f[x+h,y,z]-f[x-h,y,z])
)/(h^3)

cGradLap3D2hAnisotropic={ c6->1/2-4c1-8c2-4c3-4c4-4c5,
c11->1/2-16c1-2c10-20c2-4c3-8c4-2c5-8c7-10c8-4c9,
c12->-3+64c1+4c10+80c2+16c3+32c4+8c5+28c7+32c8+12c9 }

cGradLap3D2hIsotropic={ c5->1/6-8c1-10c2-2c3-4c4,
c6->-1/6+28c1+32c2+4c3+12c4, c8->1/32-4c1-2c2-1/4c3-2c7-1/8c10,
c9->1/48+4c1+1/2c3-2c4+2c7+1/4c10, c11->-11/48+24c1+20c2+1/2c3+8c4+4c7-7/4c10,
c12->-5/12-80c1-64c2-2c3-24c4-12c7+3c10 }

eGradLap3D2hIsotropic=1/4




(***************** 10-point O(h^2) anisotropic stencil ******************)

GradLap3D2hAniso10p[f_,x_,y_,z_,h_]:=(
(       f[x+2h,y+h,z+h]-f[x-2h,y+h,z+h]+f[x+2h,y-h,z+h]-f[x-2h,y-h,z+h]+
        f[x+2h,y+h,z-h]-f[x-2h,y+h,z-h]+f[x+2h,y-h,z-h]-f[x-2h,y-h,z-h])
-8(   f[x+h,y,z]-f[x-h,y,z])
)/(8 h^3)

cGradLap3D2hAniso10p={ c1->0, c2->0, c3->1/8, c4->0, c5->0, c6->0, 
c7->0, c8->0, c9->0, c10->0, c11->0, c12->-1 }




(***************** 12-point O(h^2) anisotropic stencil ******************)

GradLap3D2hAniso12p[f_,x_,y_,z_,h_]:=(
(	f[x+2h,y,z]-f[x-2h,y,z])
+(	f[x+h,y+h,z]-f[x-h,y+h,z]+f[x+h,y-h,z]-f[x-h,y-h,z]+
        f[x+h,y,z-h]-f[x-h,y,z-h]+f[x+h,y,z+h]-f[x-h,y,z+h])
-6(	f[x+h,y,z]-f[x-h,y,z])
)/(2h^3)

cGradLap3D2hAniso12p={ c1->0, c2->0, c3->0, c4->0, c5->0, 
c6->1/2, c7->0, c8->0, c9->0, c10->0, c11->1/2, c12->-3 }



(***************** 28-point O(h^2) isotropic stencil ********************)

GradLap3D2hIso28p[f_,x_,y_,z_,h_]:=(
(       f[x+2h,y+2h,z]-f[x-2h,y+2h,z]+f[x+2h,y-2h,z]-f[x-2h,y-2h,z]+
        f[x+2h,y,z-2h]-f[x-2h,y,z-2h]+f[x+2h,y,z+2h]-f[x-2h,y,z+2h])
+8(     f[x+2h,y,z]-f[x-2h,y,z])
+6(     f[x+h,y+h,z+h]-f[x-h,y+h,z+h]+f[x+h,y-h,z+h]-f[x-h,y-h,z+h]+
        f[x+h,y+h,z-h]-f[x-h,y+h,z-h]+f[x+h,y-h,z-h]-f[x-h,y-h,z-h])
-8(     f[x+h,y+h,z]-f[x-h,y+h,z]+f[x+h,y-h,z]-f[x-h,y-h,z]+
        f[x+h,y,z-h]-f[x-h,y,z-h]+f[x+h,y,z+h]-f[x-h,y,z+h])
-16(    f[x+h,y,z]-f[x-h,y,z])
)/(24h^3)

cGradLap3D2hIso28p={ c1->0, c2->0, c3->0, c4->1/24, c5->0, c6->1/3, 
c7->0,  c8->0, c9->0, c10->1/4, c11->-1/3, c12->-2/3 }

eGradLap3D2hIso28p=1/4



(***************** 36-point O(h^2) isotropic stencil ********************)

GradLap3D2hIso36p[f_,x_,y_,z_,h_]:=(
2(      f[x+2h,y,z+h]-f[x-2h,y,z+h]+f[x+2h,y,z-h]-f[x-2h,y,z-h]+
        f[x+2h,y+h,z]-f[x-2h,y+h,z]+f[x+2h,y-h,z]-f[x-2h,y-h,z])
-2(     f[x+2h,y,z]-f[x-2h,y,z])
+(      f[x+h,y,z+2h]-f[x-h,y,z+2h]+f[x+h,y,z-2h]-f[x-h,y,z-2h]+
        f[x+h,y+2h,z]-f[x-h,y+2h,z]+f[x+h,y-2h,z]-f[x-h,y-2h,z])
+3(     f[x+h,y+h,z+h]-f[x-h,y+h,z+h]+f[x+h,y-h,z+h]-f[x-h,y-h,z+h]+
        f[x+h,y+h,z-h]-f[x-h,y+h,z-h]+f[x+h,y-h,z-h]-f[x-h,y-h,z-h])
-8(     f[x+h,y+h,z]-f[x-h,y+h,z]+f[x+h,y-h,z]-f[x-h,y-h,z]+
        f[x+h,y,z-h]-f[x-h,y,z-h]+f[x+h,y,z+h]-f[x-h,y,z+h])
+4(     f[x+h,y,z]-f[x-h,y,z])
)/(12h^3)

cGradLap3D2hIso36p={ c1->0, c2->0, c3->0, c4->0, c5->1/6, c6->-1/6,
c7->0, c8->0, c9->1/12, c10->1/4, c11->-2/3, c12->1/3 }

eGradLap3D2hIso36p=1/4


(***************** 40-point O(h^2) isotropic stencil ********************)

GradLap3D2hIso40p[f_,x_,y_,z_,h_]:=(
16(    f[x+2h,y,z+h]-f[x-2h,y,z+h]+f[x+2h,y,z-h]-f[x-2h,y,z-h]+
        f[x+2h,y+h,z]-f[x-2h,y+h,z]+f[x+2h,y-h,z]-f[x-2h,y-h,z])
-16(    f[x+2h,y,z]-f[x-2h,y,z])
+3(	f[x+h,y+2h,z+h]-f[x-h,y+2h,z+h]+f[x+h,y-2h,z+h]-f[x-h,y-2h,z+h]+
        f[x+h,y+2h,z-h]-f[x-h,y+2h,z-h]+f[x+h,y-2h,z-h]-f[x-h,y-2h,z-h]+
	f[x+h,y+h,z+2h]-f[x-h,y+h,z+2h]+f[x+h,y-h,z+2h]-f[x-h,y-h,z+2h]+
        f[x+h,y+h,z-2h]-f[x-h,y+h,z-2h]+f[x+h,y-h,z-2h]-f[x-h,y-h,z-2h])
+2(    f[x+h,y,z+2h]-f[x-h,y,z+2h]+f[x+h,y,z-2h]-f[x-h,y,z-2h]+
        f[x+h,y+2h,z]-f[x-h,y+2h,z]+f[x+h,y-2h,z]-f[x-h,y-2h,z])
-22(   f[x+h,y+h,z]-f[x-h,y+h,z]+f[x+h,y-h,z]-f[x-h,y-h,z]+
        f[x+h,y,z-h]-f[x-h,y,z-h]+f[x+h,y,z+h]-f[x-h,y,z+h])
-40(   f[x+h,y,z]-f[x-h,y,z])
)/(96 h^3)

cGradLap3D2hIso40p={ c1->0, c2->0, c3->0, c4->0, c5->1/6, c6->-1/6, c7->0,
c8->1/32, c9->1/48, c10->0, c11->-11/48, c12->-5/12 }

eGradLap3D2hIso40p=1/4








(************************************************************************)
(*                                                                      *)
(*               CHECK ROUTINES TO CONTROL THE RESULTS                  *)
(*                                                                      *)
(************************************************************************)


(***** Check whether the explicite representations of the stencils ******)
(***** identical to the coefficient lists given *************************)

Print["Checking for differences between coefficients and explicite stencils"]
Print["Expected output: lots of zeros"]
errorList={}

errorList=Append[ errorList,
	Simplify[ ( Laplacian2D2hGeneral[f,x,y,h] /. cLaplacian2D2hAniso5p )
	- Laplacian2D2hAniso5p[f,x,y,h] ]]
errorList=Append[ errorList,
	Simplify[ ( Laplacian2D2hGeneral[f,x,y,h] /. cLaplacian2D2hIso9p )
	- Laplacian2D2hIso9p[f,x,y,h] ]]
	
errorList=Append[ errorList,
	Simplify[ ( Laplacian2D4hGeneral[f,x,y,h] /. cLaplacian2D4hAniso9p )
	- Laplacian2D4hAniso9p[f,x,y,h] ]]
errorList=Append[ errorList,
	Simplify[ ( Laplacian2D4hGeneral[f,x,y,h] /. cLaplacian2D4hIso17p )
	- Laplacian2D4hIso17p[f,x,y,h] ]]
errorList=Append[ errorList,
	Simplify[ ( Laplacian2D4hGeneral[f,x,y,h] /. cLaplacian2D4hIso21p )
	- Laplacian2D4hIso21p[f,x,y,h] ]]

errorList=Append[ errorList,
	Simplify[ ( Laplacian3D2hGeneral[f,x,y,z,h] /. cLaplacian3D2hAniso7p )
	- Laplacian3D2hAniso7p[f,x,y,z,h] ]]
errorList=Append[ errorList,
	Simplify[ ( Laplacian3D2hGeneral[f,x,y,z,h] /. cLaplacian3D2hIso15p )
	- Laplacian3D2hIso15p[f,x,y,z,h] ]]
errorList=Append[ errorList,
	Simplify[ ( Laplacian3D2hGeneral[f,x,y,z,h] /. cLaplacian3D2hIso19p )
	- Laplacian3D2hIso19p[f,x,y,z,h] ]]
errorList=Append[ errorList,
	Simplify[ ( Laplacian3D2hGeneral[f,x,y,z,h] /. cLaplacian3D2hIso21p )
	- Laplacian3D2hIso21p[f,x,y,z,h] ]]
errorList=Append[ errorList,
	Simplify[ ( Laplacian3D2hGeneral[f,x,y,z,h] /. cLaplacian3D2hIso27p )
	- Laplacian3D2hIso27p[f,x,y,z,h] ]]

errorList=Append[ errorList,
	Simplify[ ( Bilaplacian2D2hGeneral[f,x,y,h] /. cBilaplacian2D2hAniso13p)
	- Bilaplacian2D2hAniso13p[f,x,y,h] ]]
errorList=Append[ errorList,
	Simplify[ ( Bilaplacian2D2hGeneral[f,x,y,h] /. cBilaplacian2D2hIso17p )
	- Bilaplacian2D2hIso17p[f,x,y,h] ]]
errorList=Append[ errorList,
	Simplify[ ( Bilaplacian2D2hGeneral[f,x,y,h] /. cBilaplacian2D2hIso21p )
	- Bilaplacian2D2hIso21p[f,x,y,h] ]]

errorList=Append[ errorList,
	Simplify[ ( Bilaplacian2D4hGeneral[f,x,y,h] /. cBilaplacian2D4hAniso21p )
	- Bilaplacian2D4hAniso21p[f,x,y,h] ]]
errorList=Append[ errorList,
	Simplify[ ( Bilaplacian2D4hGeneral[f,x,y,h] /. cBilaplacian2D4hAniso25p )
	- Bilaplacian2D4hAniso25p[f,x,y,h] ]]
errorList=Append[ errorList,
	Simplify[ ( Bilaplacian2D4hGeneral[f,x,y,h] /. cBilaplacian2D4hIso33p )
	- Bilaplacian2D4hIso33p[f,x,y,h] ]]
errorList=Append[ errorList,
	Simplify[ ( Bilaplacian2D4hGeneral[f,x,y,h] /. cBilaplacian2D4hIso37p )
	- Bilaplacian2D4hIso37p[f,x,y,h] ]]

errorList=Append[ errorList,
	Simplify[ ( Bilaplacian3D2hGeneral[f,x,y,z,h] /.cBilaplacian3D2hAniso21p )
	- Bilaplacian3D2hAniso21p[f,x,y,z,h] ]]
errorList=Append[ errorList,
	Simplify[ ( Bilaplacian3D2hGeneral[f,x,y,z,h] /.cBilaplacian3D2hAniso25p )
	- Bilaplacian3D2hAniso25p[f,x,y,z,h] ]]
errorList=Append[ errorList,
	Simplify[ ( Bilaplacian3D2hGeneral[f,x,y,z,h] /. cBilaplacian3D2hIso41p )
	- Bilaplacian3D2hIso41p[f,x,y,z,h] ]]
errorList=Append[ errorList,
	Simplify[ ( Bilaplacian3D2hGeneral[f,x,y,z,h] /. cBilaplacian3D2hIso52p )
	- Bilaplacian3D2hIso52p[f,x,y,z,h] ]]
errorList=Append[ errorList,
	Simplify[ ( Bilaplacian3D2hGeneral[f,x,y,z,h] /. cBilaplacian3D2hIso57p )
	- Bilaplacian3D2hIso57p[f,x,y,z,h] ]]

errorList=Append[ errorList,
	Simplify[ ( GradLap2D2hGeneral[f,x,y,h] /. cGradLap2D2hAniso6p )
	- GradLap2D2hAniso6p[f,x,y,h] ]]
errorList=Append[ errorList,
	Simplify[ ( GradLap2D2hGeneral[f,x,y,h] /. cGradLap2D2hAniso8p )
	- GradLap2D2hAniso8p[f,x,y,h] ]]
errorList=Append[ errorList,
	Simplify[ ( GradLap2D2hGeneral[f,x,y,h] /. cGradLap2D2hIso12p )
	- GradLap2D2hIso12p[f,x,y,h] ]]
errorList=Append[ errorList,
	Simplify[ ( GradLap2D2hGeneral[f,x,y,h] /. cGradLap2D2hIso16p )
	- GradLap2D2hIso16p[f,x,y,h] ]]

errorList=Append[ errorList,
	Simplify[ ( GradLap2D4hGeneral[f,x,y,h] /. cGradLap2D4hAniso14p )
	- GradLap2D4hAniso14p[f,x,y,h] ]]
errorList=Append[ errorList,
	Simplify[ ( GradLap2D4hGeneral[f,x,y,h] /. cGradLap2D4hAniso18p )
	- GradLap2D4hAniso18p[f,x,y,h] ]]
errorList=Append[ errorList,
	Simplify[ ( GradLap2D4hGeneral[f,x,y,h] /. cGradLap2D4hIso26p )
	- GradLap2D4hIso26p[f,x,y,h] ]]
errorList=Append[ errorList,
	Simplify[ ( GradLap2D4hGeneral[f,x,y,h] /. cGradLap2D4hIso30p )
	- GradLap2D4hIso30p[f,x,y,h] ]]

errorList=Append[ errorList,
	Simplify[ ( GradLap3D2hGeneral[f,x,y,z,h] /. cGradLap3D2hAniso10p )
	- GradLap3D2hAniso10p[f,x,y,z,h] ]]
errorList=Append[ errorList,
	Simplify[ ( GradLap3D2hGeneral[f,x,y,z,h] /. cGradLap3D2hAniso12p )
	- GradLap3D2hAniso12p[f,x,y,z,h] ]]
errorList=Append[ errorList,
	Simplify[ ( GradLap3D2hGeneral[f,x,y,z,h] /. cGradLap3D2hIso28p )
	- GradLap3D2hIso28p[f,x,y,z,h] ]]
errorList=Append[ errorList,
	Simplify[ ( GradLap3D2hGeneral[f,x,y,z,h] /. cGradLap3D2hIso36p )
	- GradLap3D2hIso36p[f,x,y,z,h] ]]
errorList=Append[ errorList,
	Simplify[ ( GradLap3D2hGeneral[f,x,y,z,h] /. cGradLap3D2hIso40p )
	- GradLap3D2hIso40p[f,x,y,z,h] ]]

Print[errorList]
Print[]


(***** Create polynomials that we need to check whether the stencils *****)
(***** indeed compute the desired differential expressions ***************)

Clear[i,j,a2D,u2D];
n=10;
u2D=0;
Do[If[i+j<=n,u2D=u2D+a2D[i,j] x^i y^j,Continue[]],
{i,0,n},{j,0,n}];
f2D[x_,y_]=u2D;
lap2D=D[u2D,{x,2}]+D[u2D,{y,2}];
lapp2D=D[lap2D,{x,2}]+D[lap2D,{y,2}];
lappp2D=D[lapp2D,{x,2}]+D[lapp2D,{y,2}];
lapppp2D=D[lappp2D,{x,2}]+D[lappp2D,{y,2}];
glap2D=D[lap2D,x];
glapp2D=D[lapp2D,x];
glappp2D=D[lappp2D,x];

Clear[i,j,l,a3D,u3D];
n=8;
u3D=0;
Do[If[i+j+l<=n,u3D=u3D+a3D[i,j,l] x^i y^j z^l,Continue[]],
{i,0,n},{j,0,n},{l,0,n}];
f3D[x_,y_,z_]=u3D;
lap3D=D[u3D,{x,2}]+D[u3D,{y,2}]+D[u3D,{z,2}];
lapp3D=D[lap3D,{x,2}]+D[lap3D,{y,2}]+D[lap3D,{z,2}];
lappp3D=D[lapp3D,{x,2}]+D[lapp3D,{y,2}]+D[lapp3D,{z,2}];
glap3D=D[lap3D,x];
glapp3D=D[lapp3D,x];
glappp3D=D[lappp3D,x];


(***** First we check all anisotropic expressions since for those we *****)
(***** only need to check the result up to given order, not worrying *****)
(***** about the form of the error term **********************************)

Print["Check of the error order for the anisotropic stencils"]
Print["Expected Output: h^2 h^4 h^2  h^2 h^4 h^4 h^2 h^2  h^2 h^2 h^4 h^4 h^2 h^2"]
errorList={}
errorList=Append[errorList,
Series[ Laplacian2D2hAniso5p[f2D,x,y,h] - lap2D, {h,0,1}]]
errorList=Append[errorList,
Series[ Laplacian2D4hAniso9p[f2D,x,y,h] - lap2D, {h,0,3}]]
errorList=Append[errorList,
Series[ Laplacian3D2hAniso7p[f3D,x,y,z,h] - lap3D, {h,0,1}]]
errorList=Append[errorList,

Series[ Bilaplacian2D2hAniso13p[f2D,x,y,h] - lapp2D, {h,0,1}]]
errorList=Append[errorList,
Series[ Bilaplacian2D4hAniso21p[f2D,x,y,h] - lapp2D, {h,0,3}]]
errorList=Append[errorList,
Series[ Bilaplacian2D4hAniso25p[f2D,x,y,h] - lapp2D, {h,0,3}]]
errorList=Append[errorList,
Series[ Bilaplacian3D2hAniso21p[f3D,x,y,z,h] - lapp3D, {h,0,1}]]
errorList=Append[errorList,
Series[ Bilaplacian3D2hAniso25p[f3D,x,y,z,h] - lapp3D, {h,0,1}]]

errorList=Append[errorList,
Series[ GradLap2D2hAniso6p[f2D,x,y,h] - glap2D, {h,0,1}]]
errorList=Append[errorList,
Series[ GradLap2D2hAniso8p[f2D,x,y,h] - glap2D, {h,0,1}]]
errorList=Append[errorList,
Series[ GradLap2D4hAniso14p[f2D,x,y,h] - glap2D, {h,0,3}]]
errorList=Append[errorList,
Series[ GradLap2D4hAniso18p[f2D,x,y,h] - glap2D, {h,0,3}]]
errorList=Append[errorList,
Series[ GradLap3D2hAniso10p[f3D,x,y,z,h] - glap3D, {h,0,1}]]
errorList=Append[errorList,
Series[ GradLap3D2hAniso12p[f3D,x,y,z,h] - glap3D, {h,0,1}]]


Print[errorList]
Print[]




(***** Now we check whether the isotropic results have the intended ******)
(***** error in first order. *********************************************)

Print["Checking non-specified truncation error for the Laplacian"]
Print["Expected output: h^4  h^6 h^6  h^4 h^4 h^4 h^4"]
errorList={}

errorList=Append[errorList,Simplify[ Series[ Simplify[
Laplacian2D2hIso9p[f2D,x,y,h] - lap2D 
	- eLaplacian2D2hIso9p h^2 lapp2D],{h,0,3}]]]
	
errorList=Append[errorList,Simplify[ Series[ Simplify[
Laplacian2D4hIso17p[f2D,x,y,h] - lap2D 
	- eLaplacian2D4hIso17p h^4 lappp2D],{h,0,5}]]]
errorList=Append[errorList,Simplify[Series[ Simplify[
Laplacian2D4hIso21p[f2D,x,y,h] - lap2D 
	- eLaplacian2D4hIso21p h^4 lappp2D],{h,0,5}]]]

errorList=Append[errorList,Simplify[ Series[ Simplify[
Laplacian3D2hIso15p[f3D,x,y,z,h] - lap3D 
	- eLaplacian3D2hIso15p h^2 lapp3D],{h,0,3}]]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
Laplacian3D2hIso19p[f3D,x,y,z,h] - lap3D 
	- eLaplacian3D2hIso19p h^2 lapp3D],{h,0,3}]]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
Laplacian3D2hIso21p[f3D,x,y,z,h] - lap3D 
	- eLaplacian3D2hIso21p h^2 lapp3D],{h,0,3}]]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
Laplacian3D2hIso27p[f3D,x,y,z,h] - lap3D 
	- eLaplacian3D2hIso27p h^2 lapp3D],{h,0,3}]]]
	
Print[errorList]
Print[]
Print["Checking non-specified truncation error for Bilaplacian"]
Print["Expected output: h^4 h^4  h^6 h^6  h^4 h^4 h^4"]
errorList={}

errorList=Append[errorList,Simplify[ Series[ Simplify[
Bilaplacian2D2hIso17p[f2D,x,y,h] - lapp2D 
	- eBilaplacian2D2hIso17p h^2 lappp2D],{h,0,3}]]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
Bilaplacian2D2hIso21p[f2D,x,y,h] - lapp2D 
	- eBilaplacian2D2hIso21p h^2 lappp2D],{h,0,3}]]]

errorList=Append[errorList,Simplify[ Series[ Simplify[
Bilaplacian2D4hIso33p[f2D,x,y,h] - lapp2D 
	- eBilaplacian2D4hIso33p h^4 lapppp2D],{h,0,5}]]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
Bilaplacian2D4hIso37p[f2D,x,y,h] - lapp2D 
	- eBilaplacian2D4hIso37p h^4 lapppp2D],{h,0,5}]]]

errorList=Append[errorList,Simplify[ Series[ Simplify[
Bilaplacian3D2hIso41p[f3D,x,y,z,h] - lapp3D 
	- eBilaplacian3D2hIso41p h^2 lappp3D],{h,0,3}]]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
Bilaplacian3D2hIso52p[f3D,x,y,z,h] - lapp3D 
	- eBilaplacian3D2hIso52p h^2 lappp3D],{h,0,3}]]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
Bilaplacian3D2hIso57p[f3D,x,y,z,h] - lapp3D 
	- eBilaplacian3D2hIso57p h^2 lappp3D],{h,0,3}]]]

Print[errorList]
Print[]
Print["Checking non-specified truncation error for Gradient of Laplacian"]
Print["Expected output: h^4 h^4  h^6 h^6  h^4 h^4 h^4"]
errorList={}

errorList=Append[errorList,Simplify[ Series[ Simplify[
GradLap2D2hIso12p[f2D,x,y,h] - glap2D 
	- eGradLap2D2hIso12p h^2 glapp2D],{h,0,3}]]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
GradLap2D2hIso16p[f2D,x,y,h] - glap2D 
	- eGradLap2D2hIso16p h^2 glapp2D],{h,0,3}]]]

errorList=Append[errorList,Simplify[ Series[ Simplify[
GradLap2D4hIso26p[f2D,x,y,h] - glap2D 
	- eGradLap2D4hIso26p h^4 glappp2D],{h,0,5}]]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
GradLap2D4hIso30p[f2D,x,y,h] - glap2D 
	- eGradLap2D4hIso30p h^4 glappp2D],{h,0,5}]]]

errorList=Append[errorList,Simplify[ Series[ Simplify[
GradLap3D2hIso28p[f3D,x,y,z,h] - glap3D 
	- eGradLap3D2hIso28p h^2 glapp3D],{h,0,3}]]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
GradLap3D2hIso36p[f3D,x,y,z,h] - glap3D 
	- eGradLap3D2hIso36p h^2 glapp3D],{h,0,3}]]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
GradLap3D2hIso40p[f3D,x,y,z,h] - glap3D 
	- eGradLap3D2hIso40p h^2 glapp3D],{h,0,3}]]]

Print[errorList]
Print[]


(***** Now we check the general results **********************************)

Print["Checking the general results - this can take quite some time"]
Print["Expected results: h^2 h^4  h^4 h^6  h^2 h^4"]
errorList={}

errorList=Append[errorList,
Series[ ( Laplacian2D2hGeneral[f2D,x,y,h] /. cLaplacian2D2hAnisotropic )
	- lap2D, {h,0,1}]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
( Laplacian2D2hGeneral[f2D,x,y,h] /. cLaplacian2D2hIsotropic ) - lap2D 
	- eLaplacian2D2hIsotropic h^2 lapp2D],{h,0,3}]]]

errorList=Append[errorList,
Series[ ( Laplacian2D4hGeneral[f2D,x,y,h] /. cLaplacian2D4hAnisotropic )
	- lap2D, {h,0,3}]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
( Laplacian2D4hGeneral[f2D,x,y,h] /. cLaplacian2D4hIsotropic ) - lap2D 
	- eLaplacian2D4hIsotropic h^4 lappp2D],{h,0,5}]]]

errorList=Append[errorList,
Series[ ( Laplacian3D2hGeneral[f3D,x,y,z,h] /. cLaplacian3D2hAnisotropic )
	- lap3D, {h,0,1}]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
( Laplacian3D2hGeneral[f3D,x,y,z,h] /. cLaplacian3D2hIsotropic ) - lap3D 
	- eLaplacian3D2hIsotropic h^2 lapp3D],{h,0,3}]]]

Print[errorList]
errorList={}

errorList=Append[errorList,
Series[ ( Bilaplacian2D2hGeneral[f2D,x,y,h] /. cBilaplacian2D2hAnisotropic )
	- lapp2D, {h,0,1}]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
( Bilaplacian2D2hGeneral[f2D,x,y,h] /. cBilaplacian2D2hIsotropic ) - lapp2D 
	- eBilaplacian2D2hIsotropic h^2 lappp2D],{h,0,3}]]]

errorList=Append[errorList,
Series[ ( Bilaplacian2D4hGeneral[f2D,x,y,h] /. cBilaplacian2D4hAnisotropic )
	- lapp2D, {h,0,3}]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
( Bilaplacian2D4hGeneral[f2D,x,y,h] /. cBilaplacian2D4hIsotropic ) - lapp2D 
	- eBilaplacian2D4hIsotropic h^4 lapppp2D],{h,0,5}]]]

errorList=Append[errorList,
Series[ ( Bilaplacian3D2hGeneral[f3D,x,y,z,h] /. cBilaplacian3D2hAnisotropic )
	- lapp3D, {h,0,1}]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
( Bilaplacian3D2hGeneral[f3D,x,y,z,h] /. cBilaplacian3D2hIsotropic ) - lapp3D 
	- eBilaplacian3D2hIsotropic h^2 lappp3D],{h,0,3}]]]

Print[errorList]
errorList={}

errorList=Append[errorList,
Series[ ( GradLap2D2hGeneral[f2D,x,y,h] /. cGradLap2D2hAnisotropic )
	- glap2D, {h,0,1}]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
( GradLap2D2hGeneral[f2D,x,y,h] /. cGradLap2D2hIsotropic ) - glap2D 
	- eGradLap2D2hIsotropic h^2 glapp2D],{h,0,3}]]]

errorList=Append[errorList,
Series[ ( GradLap2D4hGeneral[f2D,x,y,h] /. cGradLap2D4hAnisotropic )
	- glap2D, {h,0,3}]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
( GradLap2D4hGeneral[f2D,x,y,h] /. cGradLap2D4hIsotropic ) - glap2D 
	- eGradLap2D4hIsotropic h^4 glappp2D],{h,0,5}]]]

errorList=Append[errorList,
Series[ ( GradLap3D2hGeneral[f3D,x,y,z,h] /. cGradLap3D2hAnisotropic )
	- glap3D, {h,0,1}]]
errorList=Append[errorList,Simplify[ Series[ Simplify[
( GradLap3D2hGeneral[f3D,x,y,z,h] /. cGradLap3D2hIsotropic ) - glap3D 
	- eGradLap3D2hIsotropic h^2 glapp3D],{h,0,3}]]]
