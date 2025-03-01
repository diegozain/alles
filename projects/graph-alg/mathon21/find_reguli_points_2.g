# Read("find_reguli_points.g"); 
# 
# 
pg := ProjectiveSpace(2,2);
# 
G := CollineationGroup(pg);
# 
# 
# 
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
# 
result := [P[1]];
# 
update := function(Lista)
# 
# P \ Lista
# 
  local j,i,PP;
  PP := ShallowCopy(P);
  for i in [1..Length(Lista)] do
    j := Position(PP,Lista[i]);
    Remove(PP,j);
    od;
  return PP;
end;;
# 
P_possible := update(result);

P_possible_fixed := P_possible;

find_reguli_points := function(i, size_of_reguli, alpha)
# 
# performs backtrack on P to find an SPG-reguli set
# 
  local j, what_check_was;
#   
#   point = choose_candidate(P_possible);
# 
  if Length(result) = size_of_reguli then 
    break; 
  fi;
  if check_alpha(P_possible[i], result, alpha) = false then
      what_check_was := false;
  elif check_alpha(P_possible[i], result, alpha) = true then
      what_check_was := true;
      Add(result, P_possible[i]);
      P_possible := update(result);
  fi;
  for j in [2..Length(P_possible_fixed)] do
    if what_check_was = false then
      find_reguli_points(i+1, size_of_reguli, alpha);
    elif what_check_was = true then
      find_reguli_points(i, size_of_reguli, alpha);
    fi;
  od;
  return result;
end;;

# 

check_alpha := function(point, Points, alpha)
# 
# checks if "point" satisfies alpha property with "Points"
# 
# for each line "point-Points[i]", checks if that line has more than
# "alpha-1" points of "Points"
# 
# in PG(2,2^k) there are 2^k + 2 pts in the reguli ((q+1)+1) for alpha = 1.
# in PG(2,2^k) there are 3*2^k + 4 pts in the reguli (3(q+1)+1) for alpha = 3.
# 
  local p, line, pp, i, check;
  check := true;
  for p in Points do
    i := 0;
    line := Span(point, p);
    for pp in Points do
      if IsIncident(pp, line) = true then
	i := i+1;
      fi;
    od;      
	if i>alpha then
	  check := false;
	  break;
	fi;
  od;
    return check;
end;;


# gap> find_reguli_points(1,6,1);                                                                   

# gap> find_reguli_points(1,16,3);
