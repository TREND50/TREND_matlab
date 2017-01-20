function y = Agglomerate( x, g, varargin )
% AGGLOMERATE Agglomerates a binary data set.
%   Y = AGGLOMERATE(X,G) returns the aglomerated binary
%   data set built from the set X, with granularity G.
%   Isolated '1' in X, at more than G from others are
%   deleted. Holes of width smaller than G are filled.

  x = x(:)';
  
  g1 = g;
  if ~isempty( varargin )
    g1 = varargin{ 1 };
  end
  
  k = g + 1;
  k1 = g1 + 1;
  
  % Filter isolated 1
  %===
  if g1 > 1
    xx = [ x, zeros( 1, 2*k1+1 ) ];
    xx = filter( ones( 1, 2*k1 + 1 ), 1, xx );
    x = ( xx( k1+1:end-k1-1 ) > 1 ) & ( x == 1 );
  end
    
  % Fill holes
  %===
  mask = 2.^[ k-1:-1:0 ];
  
  yb = filter( mask, 1, x ); % left side
  K  = yb  > 0;
  J  = yb == 0;
  yb( K ) = k - 1 - floor( log( yb( K ) )/log( 2 ) );
  yb( J ) = k + 1;
  
  ya = filter( mask, 1, x( end:-1:1 ) ); % right side
  ya = ya( end:-1:1 );
  K  = ya  > 0;
  J  = ya == 0;
  ya( K ) = k - 1 - floor( log( ya( K ) )/log( 2 ) );
  ya( J ) = k + 1;
  
  y = x | ( ya + yb <= k );