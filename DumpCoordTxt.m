function DumpCoordTxt( Struct )


  SharedGlobals;

  Det = Struct.Setup.Det;
  
  kfirst = min( [ Det.Name ] );
  klast  = max( [ Det.Name ] );
  
  filename=[TEXT_PATH 'coord_antennas.txt']
  fid = fopen(filename, 'w+' );
  for k = kfirst:klast
    j = find( [ Det.Name ] == k );
    if isempty( j )
      fprintf( fid, '%3d %6.1f %6.1f %6.1f\n', ...
        k, 0, 0, 0 );
    else
      fprintf( fid, '%3d %6.1f %6.1f %6.1f\n', ...
        k, Det( j ).X, Det( j ).Y, Det( j ).Z );
    end
  end
  fclose( fid );

 filename=[TEXT_PATH 'coord_detectors.txt'];
 fid = fopen(filename, 'w+' );
   for k = kfirst:klast
     j = find( [ Det.Name ] == k );
     if isempty( j )
       fprintf( fid, '%3d %6.1f %6.1f %6.1f\n', ...
         k, 0, 0, 0 );
     else
       fprintf( fid, '%3d %6.1f %6.1f %6.1f\n', ...
       k, Det( j ).X, Det( j ).Y, Det( j ).Z );
     end
 end
 fclose( fid );

