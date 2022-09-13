% This version assumes that nodes are numbered 1,2,...,n
function N = loadnetwork(filename_nod, filename_edg)
  A = load(filename_nod);
  N.n = size(A, 1);
  N.node_idx = A(:,1);
  N.P = A(:,2:4);
  N.valency = zeros(size(N.node_idx));

  N.cons = zeros(N.n, max(N.valency));
  B = load(filename_edg);
  for k=1:size(B,1)
    i = B(k, 2);
    j = B(k, 3);
    N.valency(i) = N.valency(i)+1;
    N.valency(j) = N.valency(j)+1;
    N.cons(i, N.valency(i)) = j;
    N.cons(j, N.valency(j)) = i;
  end;
end
