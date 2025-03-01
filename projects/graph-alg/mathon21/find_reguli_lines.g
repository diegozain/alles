# Read("find_reguli_lines.g");
# 
# 
pg := ProjectiveSpace(3,2);
# 
G := CollineationGroup(pg);
# 
# 
# 
# 
L_0 := AsList(Lines(pg));
L := [Length(L_0)];
for k in [1..Length(L_0)] do
  L[k] := L_0[k];
od;
# 
P_0 := AsList(Points(pg));
P := [Length(P_0)];
for k in [1..Length(P_0)] do
  P[k] := P_0[k];
od;
# 
# 
# 
# 
# 
result := [L[1]];
# 
LINES := ShallowCopy(L);
# 
IsIncident_mio := function(ele, elem)
local meet;
meet := Meet(ele, elem);
if ProjectiveDimension(meet) = -1 then
  return false;
  else
  return true;
fi;
end;;
# 
# counter_update := 0;
# 
update_0 := function(ll,line_in_Lista)
#
# L \ [lines incident to line_in_lista]
#
  local inci;
  if IsIncident_mio(LINES[ll],line_in_Lista) = false then
    inci := false;
  elif IsIncident_mio(LINES[ll],line_in_Lista) = true then
    inci := true;
    Remove(LINES,ll);
    counter_update := counter_update + 1;
  fi;
  while counter_update<Length(LINES) do
    if inci = false then
      update_0(ll+1,line_in_Lista);
    elif inci = true then
      update_0(ll,line_in_Lista);
    fi;    
  od;
  return LINES;
end;;
#
update := function(Lista)
#
# L \ [lines incident to lines in Lista]
#
local l;
  for l in Lista do
    counter_update := 0;
    update_0(1,l);
  od;
end;;
#
check_alpha := function(Lista, alpha)
#
# checks alpha condition on Lista
#
local l,ll,plane,counter_alpha,check;
  counter_alpha := 0;
  check := true;
  for l in Lista do
    for plane in Planes_of_line(l) do
      for ll in Lista do
	if meet_mio(ll,plane_of_line(l))=0 then
	  counter_alpha = counter_alpha + 1;
	  elif meet_mio(ll,plane_of_line(l))>0 then
	    check := false;
	  break;
	fi;
      od;
    if counter_alpha > alpha then
      check := false;
      break;
    fi;
    od;
  od;
return check;
end;;
#
Planes_of_line := function(line)
#
# gives all planes through line
#
local point_list, p, Planes,plane;
  point_list := Point_list_for(line);
  for p in point_list do
    plane := Span(line, p);
    Add(Planes,plane);
  od;
  return Planes;
end;;
#  
# Some_point_list := P \ [points in lines of reguli]
PP := ShallowCopy(P);
for c in [1..Length(P)] do
  for x in L do
    if IsIncident(P[c],x) = true then
      Remove(PP,c);
    fi;
  od;
od;
Some_point_list := ShallowCopy(PP);
# 
Point_list_for := function(line)
# 
# list of points such that <point, line> is a different plane
# for every point in the list
# 
local p, plane, Point_list, i;
  p := Some_point_list[1];
  plane := Span(p, line);
  Add(Point_list, p);
  Remove(Some_point_list, 1);
  for i in [1..Length(Some_point_list)] do
    if IsIncident(Some_point_list[i], plane) = true then
      Remove(Some_point_list, i);
    fi;
  od;
  while Length(Some_point_list)>0 do
    Point_list_for(line);
  od;
  return Point_list;
end;;
#
# 
# # gap> find_reguli_points(1,6,1);                                                                   
# 
# # gap> find_reguli_points(1,16,3);
