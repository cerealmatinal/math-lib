open Lacaml.D
open Plplot

let my_function x = sin x +. exp x -. log x

let matrix_a = Mat.of_array [|[|1.; 2.; 3.|]; [|4.; 5.; 6.|]; [|7.; 8.; 9.|]|]
let vector_b = Vec.of_array [|1.; 2.; 3.|]

let matrix_mult = Mat.dot matrix_a matrix_a
let vector_sum = Vec.add vector_b vector_b

let plot_graph () =
  let x = Array.init 100 (fun i -> float i *. 0.1) in
  let y = Array.map my_function x in
  let pl = plinit () in
  plenv (-10.) 10. (-10.) 10. 0 0;
  pllab "x" "y" "My Function";
  plline x y;
  plend ()

let integrate f a b n =
  let dx = (b -. a) /. float n in
  let rec loop i acc =
    if i > n then acc *. dx
    else loop (i + 1) (acc +. f (a +. float i *. dx))
  in
  loop 0 0.

let differentiate f x =
  let h = 0.00001 in
  let slope = (f (x +. h) -. f (x -. h)) /. (2. *. h) in
  slope

let interpolate x y =
  let len = Array.length x in
  let c = Array.make_matrix len len 0. in
  for i = 0 to len - 1 do
    c.(i).(0) <- y.(i)
  done;
  for j = 1 to len - 1 do
    for i = j to len - 1 do
      c.(i).(j) <- (c.(i).(j - 1) -. c.(i - 1).(j - 1)) /.
                    (x.(i) -. x.(i - j))
    done
  done;
  fun x ->
    let sum = ref c.(len - 1).(0) in
    let prod = ref 1. in
    for i = 1 to len - 1 do
      prod := !prod *. (x -. x.(len - 1 - i));
      sum := !sum +. !prod *. c.(len - 1).(i)
    done;
    !sum

let newton_raphson f df x0 tol =
  let rec loop x =
    let delta_x = -. f x /. df x in
    if abs_float delta_x < tol then x
    else loop (x +. delta_x)
  in
  loop x0

let read_csv filename =
  let ic = open_in filename in
  let rec loop acc =
    try
      let line = input_line ic in
      loop (line :: acc)
    with
      End_of_file -> List.rev acc
  in
  let lines = loop [] in
  let data = List.map (fun line -> List.map float_of_string (String.split_on_char ',' line)) lines in
  close_in ic;
  data

let mean xs = List.fold_left (+.) 0. xs /. float_of_int (List.length xs)

let variance xs =
  let mu = mean xs in
  let square x = x *. x in
  mean (List.map (fun x -> square (x -. mu)) xs)

let stdev xs = sqrt (variance xs)

let correlation xs ys =
  let n = List.length xs in
  let mx = mean xs in
  let my = mean ys in
  let sx = stdev xs in
  let sy = stdev ys in
  let numerator = ref 0. in
  for i = 0 to n - 1 do
    numerator := !numerator +. ((List.nth xs i -. mx) *. (List.nth ys i -. my))
  done;
  !numerator /. (float_of_int (n - 1) *. sx *. sy)

let read_csv filename =
  let ic = open_in filename in
  let rec read_lines acc =
    try
      let line = input_line ic in
      let values = List.map float_of_string (Str.split (Str.regexp ",") line) in
      read_lines (values :: acc)
    with End_of_file -> close_in ic; List.rev acc
  in
  read_lines []

let linear_regression xs ys =
  let n = List.length xs in
  let x_sum = List.fold_left (+.) 0. xs in
  let y_sum = List.fold_left (+.) 0. ys in
  let xy_sum = List.fold_left (+.) 0. (List.map2 ( *. ) xs ys) in
  let x_sq_sum = List.fold_left (+.) 0. (List.map (fun x -> x *. x) xs) in
  let y_sq_sum = List.fold_left (+.) 0. (List.map (fun y -> y *. y) ys) in
  let b = (float_of_int n *. xy_sum -. x_sum *. y_sum) /. (float_of_int n *. x_sq_sum -. x_sum *. x_sum) in
  let a = (y_sum -. b *. x_sum) /. float_of_int n in)
  (a, b)

let main () =
  let window = GWindow.window ~title:"Análise de dados" ~width:800 ~height:600 () in
  window
  let vbox = GPack.vbox ~packing:window
  let filename_entry = GEdit.entry ~text:"data.csv" ~width:200 ~max_length:100 ~packing:vbox
  let run_button = GButton.button ~label:"Executar" ~packing:vbox#
  let drawing_area = GMisc.drawing_area ~packing:vbox
  let surface = Cairo.Image.create Cairo.Image.ARGB32 800 600 in
  let draw () =
    let ctx = Cairo.create surface in
    Cairo.set_source_rgb ctx 1. 1. 1.;
    Cairo.paint ctx;
    Cairo.set_line_width ctx 2.;
    Cairo.set_source_rgb ctx 0. 0. 0.;
    let data = read_csv filename_entry
    let xs = List.map List.hd data in
  let ys = List.map (fun row -> List.nth row y_index) rows in
  let xs = List.map (fun row -> List.filteri (fun i _ -> i <> y_index) row) rows in
  let xt = List.map List.rev (transpose xs) in
  let n = float_of_int (List.length rows) in
  let xtx = mat_of_list (List.map (fun row -> List.map float_of_string row) (List.map (fun row -> List.map string_of_float row) xt)) in
  let xty = mat_of_list (List.map (fun row -> [float_of_string (List.nth row y_index)]) rows) in
  let beta = gemv (gesv xtx) xty in
  let predicted = List.map (fun row -> dot_product (mat_of_list [row]) beta) xs in
  let residuals = List.map2 (-.) ys predicted in
  let mse = List.fold_left (fun acc r -> acc +. r *. r) 0. residuals /. n in
  (beta, mse)