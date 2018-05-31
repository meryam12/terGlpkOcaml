(* Aci ́erie produit bandes et rouleaux m ́etalliques pendant 40
   heures par semaine au total
   Vitesses de production : 200 bandes par heure et 140 rouleaux
   par heure.
   bandes vendues 25 euros l’unit ́e et rouleaux 30 euros l’unit ́e
   March ́e limit ́e : impossible de vendre de 6000 bandes et 4000
   rouleaux 


   Modélisation
   AMPL
   Param`etres
   Vers la R ́esolution
   Modélisation d’un problème de programmation linéaire
   Ecriture dans le langage math ́ematique, ”Mise en inéquations et
   équations”
   Variables
   Param`etres
   Contraintes
   Objectif(s) 
   x1 :le nombre de bandes `a produire
   x2 :le nombre de rouleaux `a produire

   25x1+ 30x2=z[max]
   x1≥0
   x2≥0
   x1≤6000
   x2≤4000
   1/200.x1+ 1/140.x2 ≤ 40

   langage AMPL 
   set of products
   set PROD;
   # parameters
   param heures_ouvrees >= 0;
   param vitesse_production {PROD} >= 0;
   param prix_vente {PROD} >= 0;
   param vente_max {PROD} >= 0;
   # variables
   var qte_produite {p in PROD} >= 0, <= vente_max [p];
   # objective
   maximize profit :
   sum {p in PROD} qte_produite [p] * prix_vente [p];
   # constraints
   subject to production_limitee :
   sum {p in PROD}
   (qte_produite [p] / vitesse_production [p]) <= heures_ouvree



   let x1 : int (*vitesse de production*)
   let x2 : int (*prix de vente *)
   val qt_max : int 
   val heure : int


   type produit = {
   vente_max : float; 
   prix_vente : float;
   vitesse_production : float; 
   }

   (*comparer : verifie ke toutes les qte sont inferieurs a la qte maximal*)
   let comparer = if produit < qt_max then true 
   else false

   let create (nombre_d_heures : int) (produits : produit list) =
   let module LP in 
   {  variables = produits;
     variable_kinds = foreach produits (fun produit -> (produit, Bounded (0., produit.vente_max))); 
     objective = Maximise (sum produits (fun produit -> produit.prix_vente * produit)); 
     constraints = 
       [ sum produits (fun produit -> produit / produit.vitesse_production) <= nombre_d_heures ]
   }  




   (*fonction objective 
   let maxi z = (sum prd.qt) * sum (prd.prix_unit) *)






   (* Pour savoir quel direction choisir pour l'optimisation*)
   type option_direction = Minimize | Maximize

   let sum_produits = List.fold_left ( + ) 0;;
   (*val sum_list : int list -> int = <fun>*)

   type variable_kinds =
   | Continuous_variables (** continuous variable *)
   | Int_variables (** integer variable *)


   let rec multipl : Float * produit list -> produit list =
   function (_,[]) -> [0]   
         |  (x ,[produit]) -> [x * produit.prix_vente]    
         | (x,produit::xs) ->  x * produit.prix_vente :: multipl (x,xs);; 


   let rec foreach : produit list -> float =
   function [] -> 0
         | [produit] -> produit.prix_vente
         | produit::xs -> produit.prix_vente :: foreach (xs)

   cuttingStockLP :: Problem -> Set Pattern -> LP.LP Pattern Double
   cuttingStockLP problem patterns =
   LP.LP
   { LP.direction = LP.Min
   , LP.objective = sumOver patterns $ \ p -> patternCost p `dot` p
   , LP.constraints = Map.elems $ forEach (finals problem) requestConstraint
   , LP.varBounds = forEach patterns (\p -> LP.LBound 0 `dot` p)
   , LP.varTypes = forEach patterns (\p -> LP.ContVar `dot` p)
   }
   where
   requestConstraint final =
     let qty = fromIntegral $ requestedQuantity problem final in 
       (sumOver patterns $ \ pattern ->
         finalsContributedByPattern final pattern `dot` pattern
       ) `equals` qty
       $ (showFinal final)



   (*
   let create_lp linear_program dir obj constr vb vt = 
   let linear_program in 
     of_direction linear_program dir ;
     for i = 0 to ( Array.length obj ) -1  do

     done;
 *)

   let set_direction lp dir = lp.direction = dir

   let set_objective lp obj = 
   List.map term (obj,'x'i)




   let create_lp linear_program dir obj constr vb vt = 
   let linear_program in 
   set_direction linear_program dir ;
   set_objective linear_program obj ;



   let term (coef, var) = singleton coef var ;

    let termVar (c,v) = v;

      let sumOver indices fct = let add (List.map term indices) empty ;(* forme une liste de coef*var avec tous les éléments (coef,var) de la liste indices)*)

        let forEach indices f = union f (List.map termVar indices) ;

*)



type direction = Minimize | Maximize


type variable_kind =
  | Continuous (** continuous variable *)
  | Integral (** integer variable *)


type bound =
  | Free
  | Upper of float
  | Lower of float
  | Boxed of float * float 

type 'variable linear_form =
  ('variable * float) list 

type 'variable variable_description =
  { variable : 'variable;
    identifier : string;
    bounds : bound;
    kind : variable_kind;
    objective_value : float;
  }

type constraint_kind =
  | LessThan
  | MoreThan
  | EqualsTo

type 'variable constraint_description =
  { name : string;
    constraint_kind : constraint_kind;
    right_hand_side : float;
    left_hand_side : 'variable linear_form;
  }


type 'variable linear_program = {
  direction : direction ;
  objective : 'variable linear_form; 
  variables : 'variable variable_description list;
  constraints : 'variable constraint_description list;
}


let my_lp = 
  { direction = Minimize;
    objective = [(1, 2.0); (2, 3.0)];
    variables = 
      [ { variable = 1; 
          identifier = "x1";
          bounds = Free;
          kind = Continuous;
          objective_value = 32.
        }
      ];
    constraints = []
  }


type 'variable solution =
  { value : float;
    primal : (string * float) list;
    dual : (string * float) list 
  }

type 'variable outcome =  
  | Unbounded 
  | Infeasible
  | Solution of 'variable solution 


(*primal *)
open Glp
(*code du prof
  let primal =
  let prob = Glp.create_problem () in 
  let _ = Glp.add_cols prob 3 in 
  ignore (Glp.add_rows prob 4);
  for i = 1 to 3 do 
    Glp.set_col_name prob i (Printf.sprintf "x%d" i);
    Glp.set_col_bnds prob i (Glp.Types.Double (0., 10.);
  done;

  let params = Ctypes.allocate Glp.Types.Smcp.t in
  Glp.init_smcp params;  
  Glp.simplex lp params;
  for i = 0 to 3 do
    Printf.printf "x%d = %f\n" i (Glp.get_col_prim prob i)

  done;  

*)





let essaie lp =
  (* création de problème*)
  let lp = Glp.create_problem () in
  Glp.set_problem_name lp "Essaie";
  Glp.set_objective_direction lp Maximize;
  (* remplissage avec les valeurs*)
  Glp.add_rows lp 2;
  Glp.set_row_bnds lp 1 Upper 0.0 1.0;
  Glp.set_row_bnds lp 2 Upper 0.0 2.0;
  Glp.add_cols lp 2 ;
  Glp.set_col_name lp 1 "x1";
  Glp.set_col_bnds lp 1 Lower 0.0 0.0;
  Glp.set_obj_coef lp 1 0.6;
  Glp.set_col_name lp 2 "x2";
  Glp.set_col_bnds lp 2 Lower 0.0 0.0;
  Glp.set_obj_coef lp 2 0.5;
  let i.(1) = 1;
    j.(1) = 1;
    r.(1) = 1.0;
    i.(2) = 1;
    j.(2) = 2;
    r.(2) = 2.0;
    i.(3) = 2;
    j.(3) = 1;
    r.(3) = 3.0;
    i.(4) = 2;
    j.(4) = 2;
    r.(4) = 1.0;
  in 
  Glp.load_matrix lp i j r;
  (*Solve*)
  let params = Ctypes.allocate Glp.Types.Smcp.t in
  Glp.init_smcp params;  
  Glp.simplex lp params;
  let z = Glp.get_obj_val lp; 
    let x1 = Glp.get_col_prim(lp, 1);
      let x2 = Glp.get_col_prim(lp, 2);
        Printf.printf "Z = %f\n" z
          Printf.printf "x1 = %f\n" x1
          Printf.printf "x2 = %f\n" x2
(* résolution du problème 
   maximize 0.6x1 + 0.5 x2
   subject to x1 + 2x2 <= 1
          3x1 + x2 <= 2 
*)


































