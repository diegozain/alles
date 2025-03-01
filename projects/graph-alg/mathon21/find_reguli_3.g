# Read("find_reguli_3.g");
# 
# 
pg := ProjectiveSpace(2,8);
# 
G := CollineationGroup(pg);
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
POINTS := ShallowCopy(P);
# 
counter := 0;
# 
result := [P[1]];
# 
find_reguli_points := function(i, size_reguli, alpha)
# 
# performs backtrack on P to find an SPG-reguli set
# 
  local alpha_check;
#   
#   point = choose_candidate(P_possible);
# 
  if counter>Length(POINTS) then
    POINTS := ShallowCopy(P);
    result := [P[1]];
    counter := counter+1;
    i := counter;    
  fi;  
  if i>Length(POINTS) then
    POINTS := ShallowCopy(P);
    result := [P[1]];
    counter := counter+1;
    i := counter;
  fi;
  if check_alpha(POINTS[i], result, alpha) = false then
    alpha_check := false;
  elif check_alpha(POINTS[i], result, alpha) = true then
    alpha_check := true;
    Add(result,POINTS[i]);
    Remove(POINTS,i);
#     counter := counter+1;
  fi;
  
  while Length(result)<size_reguli do
    if alpha_check = false then
      find_reguli_points(i+1, size_reguli, alpha);
    elif alpha_check = true then
      find_reguli_points(i, size_reguli, alpha);
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
    if ProjectiveDimension(line) = 0 then 
      check := false;
      break; 
    fi;
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
