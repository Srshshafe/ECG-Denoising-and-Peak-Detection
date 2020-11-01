function peaks = PeakDetection(x,ff,varargin)
N = length(x);
peaks = zeros(1,N);
thresh = .7;
rng = floor(thresh/ff);
if(nargin==3),
    flag = varargin{1};
else
    flag = abs(max(x))>abs(min(x));
end
if(flag)
    for j = 1:N,
        if(j>rng && j<N-rng)
            index = j-rng:j+rng;
        elseif(j>rng)
            index = N-2*rng:N;
        else
            index = 1:2*rng;
        end
        if(max(x(index))==x(j))
            peaks(j) = 1;
        end
    end
else
    for j = 1:N,
        if(j>rng && j<N-rng)
            index = j-rng:j+rng;
        elseif(j>rng)
            index = N-2*rng:N;
        else
            index = 1:2*rng;
        end
        if(min(x(index))==x(j))
            peaks(j) = 1;
        end
    end
end
I = find(peaks);
d = diff(I);
peaks(I(d<rng))=0;