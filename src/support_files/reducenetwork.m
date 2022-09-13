%
% N0         Network that connects ALL neighboring nodes.
% valency    Vector of target valency for all nodes.
% nAttempts  Number of attempts to modify the network
%
function [N, Es] = reducenetwork(N0, valency, nAttempts)
  if length(valency) ~= N0.n
    error('Nr of target valencies must be equal to nr of network nodes.');
  end
  if any(valency > N0.valency)
    error('No target valency can be greater than the number of neighbors.');
  end

  % Copy structure
  N = N0;

  % Remove all connections, unless they are same as target valency
  N.valency = zeros(size(N0.valency));
  N.cons = zeros(size(N0.cons));
  nonzeros = find(N0.valency > 0);

  % Calculate square difference between valencies and target valencies
  E0   = sum((N.valency - valency).^2);
  E    = E0;
  T    = 10*E0/N.n;
  sf_T = 1+log(10/N.n)/nAttempts;
  Es   = zeros(nAttempts,1);

  for k=1:nAttempts
    Es(k) = E;
    T = T*sf_T;

    % Random connection in fully populated graph
    i = nonzeros(randi(length(nonzeros))); % from node
    j = N0.cons(i,randi(N0.valency(i)));   % to node

    % Current residuals of connected nodes
    res_i = N.valency(i) - valency(i);
    res_j = N.valency(j) - valency(j);

    % Toggle the connection
    areConnected = any(N.cons(i,1:N.valency(i)) == j);
    if areConnected
      res_i_mod = res_i - 1;
      res_j_mod = res_j - 1;
    else % Nodes are not connected
      res_i_mod = res_i + 1;
      res_j_mod = res_j + 1;
    end
    dE = res_i_mod^2 + res_j_mod^2 - res_i^2 - res_j^2;
    accept = rand < exp(-dE/T);
    if accept
      if areConnected
        % Remove connection from i to j.
        j0 = find(N.cons(i,1:N.valency(i)) == j);
        N.cons(i,j0) = N.cons(i,N.valency(i));
        N.cons(i,N.valency(i)) = 0;
        N.valency(i) = N.valency(i) - 1;
        i0 = find(N.cons(j,1:N.valency(j)) == i);
        N.cons(j,i0) = N.cons(j,N.valency(j));
        N.cons(j,N.valency(j)) = 0;
        N.valency(j) = N.valency(j) - 1;
      else
        % Add connection from i to j.
        N.valency(i) = N.valency(i) + 1;
        N.cons(i, N.valency(i)) = j;
        N.valency(j) = N.valency(j) + 1;
        N.cons(j, N.valency(j)) = i;
      end
      E = E + dE;
    end
  end % the loop
end

