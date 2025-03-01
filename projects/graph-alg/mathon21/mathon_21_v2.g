#Read("mathon_21_v2.g");
ps := ProjectiveSpace(5,3);
# 
L1 := [[0,1,0,0,0,0],\
       [1,0,0,0,0,0]]\
      *Z(3)^0;
L2 := [[1,0,1,0,0,0],\
       [0,1,0,0,1,0]]\
      *Z(3)^0;
L3 := [[2,0,1,0,0,0],\
       [0,2,0,0,1,0]]\
      *Z(3)^0;
L4 := [[1,0,0,1,0,0],\
       [0,2,0,0,0,1]]\
      *Z(3)^0;
L5 := [[2,0,0,1,0,0],\
       [0,1,0,0,0,1]]\
      *Z(3)^0;
L6 := [[1,1,1,1,0,0],\
       [2,1,0,0,2,1]]\
      *Z(3)^0;
L7 := [[1,2,1,1,0,0],\
       [0,1,1,0,2,1]]\
      *Z(3)^0;
L8 := [[1,1,2,1,0,0],\
       [0,1,1,0,1,1]]\
      *Z(3)^0;
L9 := [[1,2,2,1,0,0],\
       [1,1,0,0,1,1]]\
      *Z(3)^0;
L10 := [[2,1,0,1,1,0],\
        [2,2,1,0,0,1]]\
      *Z(3)^0;
L11 := [[2,2,0,1,1,0],\
        [1,2,1,0,0,1]]\
      *Z(3)^0;
L12 := [[1,0,1,1,1,0],\
        [2,1,1,1,0,1]]\
      *Z(3)^0;
L13 := [[2,2,1,1,1,0],\
        [1,0,1,2,0,1]]\
      *Z(3)^0;
L14 := [[1,0,2,1,1,0],\
        [1,1,1,2,0,1]]\
      *Z(3)^0;
L15 := [[2,1,2,1,1,0],\
        [2,0,1,1,0,1]]\
      *Z(3)^0;
L16 := [[2,1,0,2,1,0],\
        [1,0,2,2,0,1]]\
      *Z(3)^0;
L17 := [[2,2,0,2,1,0],\
        [2,0,2,1,0,1]]\
      *Z(3)^0;
L18 := [[1,0,1,2,1,0],\
        [2,2,2,2,0,1]]\
      *Z(3)^0;
L19 := [[1,1,1,2,1,0],\
        [2,1,2,0,0,1]]\
      *Z(3)^0;
L20 := [[1,0,2,2,1,0],\
        [1,2,2,1,0,1]]\
      *Z(3)^0;
L21 := [[1,2,2,2,1,0],\
        [1,1,2,0,0,1]]\
      *Z(3)^0;
# 
21_lines := [L1, L2, L3, L4, L5, L6, L7,\
	     L8, L9, L10, L11, L12, L13, L14,\
	     L15, L16, L17, L18, L19, L20, L21];
# assign lines to PG(5,3)
line := [];
for i in [1 .. Length(21_lines)] do
  line[i] := \
  VectorSpaceToElement( ps,21_lines[i] );
od;
# 
line_list := function()
local line_listt;
line_listt := [];
  for i in [1 .. Length(21_lines)] do
    line_listt[i] := \
    VectorSpaceToElement( ps,21_lines[i] );
  od;
  return line_listt;
end;;
# 
# 7 solids
# first, all solids in which one line L (line[L]), can be.
L_solids := function(L)
local solid_list, total,largo,solid, intersec,j;
      solid_list := [];
      total := [];
      largo := Length(line);
  for j in [1 .. largo] do
    solid_list := [];
    if j = L then continue;
    else
      solid := Span(line[L],line[j]);
      for i in [1 .. largo] do
	intersec := Meet(solid, line[i]);
	if ProjectiveDimension(intersec) = 1 then
# 	  Print(" ",i," ");
	  Add(solid_list, i);
	fi;
      od;
    fi;
#     Print("\n");
    Add(total, solid_list);
  od;
  return total;
end;;
# 
loop_1 := function(L,j)
local new_line, intersec, solid, k,new_list;
      new_line := line_list();
      new_list := [];
      solid := Span(new_line[L],new_line[j]);
      k := Length(new_line);
    for i in [1..k] do
      intersec := Meet(solid, new_line[i]);
      if i = L then continue;
	elif i = j then continue;
	elif ProjectiveDimension(intersec) = 1 then
	  Add(new_list, L);
	  Add(new_list, j);
	  Add(new_list, i);
# 	update list!
      fi;
    od;
    return new_list;
end;;
# 
loop_2 := function(L)
local new_line,k,j, new_list;
new_line := line_list();
new_list := [];
k := Length(new_line);
    for j in [1..k] do
      if Length(loop_1(L,j))>0 then
	Add(new_list, loop_1(L,j));
      fi;
#       Append(new_list,loop_1(L,j));
    od;
    return new_list;
end;;
#
loop_3 := function()
local L;
  for L in [1..Length(line)] do
    Print(L, " ", loop_2(L),"\n","\n");
  od;
end;;
# 
imprime_solidos := function(lista)
local i;
  for i in lista do
    Print("L",i,"\n");
    Display(line[i]);
  od;
end;;
# 
imprime_seven_solidos := function(lista)
local i,j;
j := 1;
  for i in lista do
    Print("S",j,"\n");
    imprime_solidos(i);
    j := j+1;
  od;
end;;
7_solids := [ [ 1, 4 ], [ 2, 13 ], [ 3, 16 ],\
	      [ 6, 8 ], [ 7, 11 ], [ 9, 12 ],\
	      [ 10, 14 ] ];
solids := [];
for a in 7_solids do
  Add(solids, Span(line[a[1]],line[a[2]]));
od;
print_seven_S := function()
local i,j;
j := 1;
  for i in solids do
    Print("S",j,"\n");
    Display(i);
    j := j+1;
  od;
end;;
L_I := [[1,1,0,0,0,1],\
       [0,0,0,1,0,2]]\
      *Z(3)^0;
line_intersec := VectorSpaceToElement( ps,L_I );
# PLANES
P1  := [[1,0,0,0,0,0],\
       [0,1,0,0,0,0],\
       [0,0,0,1,0,0]]
      *Z(3)^0;
P2  := [[1,0,2,0,0,1],\
       [0,1,1,0,0,0],\
       [0,0,0,1,0,2]]
      *Z(3)^0;
P3  := [[1,0,1,0,0,1],\
       [0,1,2,0,0,0],\
       [0,0,0,1,0,2]]
      *Z(3)^0;
# assign planes to PG(5,3)
plane := [];
3_planes := [P1,P2,P3];
for i in [1 .. Length(3_planes)] do
  plane[i] := \
  VectorSpaceToElement( ps,3_planes[i] );
od;
# 
# 
#
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#Read("mathon_21_v2.g");