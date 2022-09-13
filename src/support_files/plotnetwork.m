% This version assumes that nodes are numbered 1,2,...,n
function plotnetwork(N)
  for i=1:N.n
    for nidx=1:N.valency(i)
      j = N.cons(i, nidx);
      if (i < j)
	plot3(N.P([i,j],1),N.P([i,j],2),N.P([i,j],3),'b');
        hold on;
      end
    end
  end
  hold off;
end
