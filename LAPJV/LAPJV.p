function LAPJV (n: integer; c: mat; var x,y,u,v: vec): integer;
{ as published in
  R. Jonker and A. Volgenant, University of Amsterdam,
  A Shortest Augmenting Path Algorithm
  for Dense and Sparse Linear Assignment Problems,
  Computing 38, 325-340 (1987). }
{ n: problem size;
  c: costs matrix;
  x: columns assigned to rows ;
  y: rows assigned to columns ;
  u: dual row variables ;
  v: dual column variables }

label augment;
const inf=1000000; { inf is a suitably large number }
var   f,h,i,j,k,f0,i1,j1,j2,u1,u2,min,last,low,up: integer;
      col,d,free,pred: vec;

{ col  : array of columns, scanned                 (k=1..low-1),
                           labeled and unscanned   (k=low..up-1),
                           unlabeled                (k=up..n);
  d    : shortest path lengths;
  free : unassigned rows (number f0, index f);
  pred : predecessor-array for shortest path tree;
  i,i1 : row indices;  j,j1,j2: column indices;
  last : last column in col-array with d[j]<min. }

begin
   for i:=1 to n do x[i]:=0;

   for j:=n downto 1 do  { #### COLUMN REDUCTION }
   begin
      col[j]:=j; h:=c[1,j]; i1:=1;
      for i:=2 to n do if c[i,j]<h then begin h:=c[i,j]; i1:=i end;
      v[j]:=h;
      if x[i1]=0 then begin x[i1]:=j; y[j]:=i1 end
         else begin x[i1]:=-abs(x[i1]); y[j]:=0 end
   end;

   f:=0; { #### REDUCTION TRANSFER }
   for i:=1 to n do
   if x[i]=0 then       { ## unassigned row in free-array }
      begin f:=f+1; free[f]:=i end else
   if x[i]<0 then       { ## no reduction transfer possible }
      x[i]:=-x[i] else  { ## reduction transfer from assigned row }
   begin 
      j1:=x[i]; min:=inf; 
      for j:=1 to n do if j<>j1 then 
         if c[i,j]-v[j]<min then min:=c[i,j]-v[j]; 
      v[j1]:=v[j1]-min 
   end; 
 
   cnt:=0;{ #### AUGMENTING ROW REDUCTION } 
   repeat 
      k:=1; f0:=f; f:=0; 
      while k<=f0 do
      begin 
         i:=free[k]; k:=k+1; u1:=c[i,1]-v[1]; j1:=1; u2:=inf; 
         for j:=2 to n do 
         begin 
            h:=c[i,j]-v[j]; 
            if h<u2 then 
               if h>=u1 then begin u2:=h; j2:=j end 
                  else begin u2:=u1; u1:=h; j2:=j1; j1:=j end 
         end; 
         i1:=y[j1]; 
         if u1<u2 then v[j1]:=v[j1]-u2+u1
            else if i1>0 then begin j1:=j2; i1:=y[j1] end; 
         if i1>0 then 
            if u1<u2 then begin k:=k-1; free[k]:=i1 end 
               else begin f:=f+1; free[f]:=i1 end; 
         x[i]:=j1; y[j1]:=i 
      end; 
      cnt:=cnt+1 
   until cnt=2; { ## routine applied twice } 
 
   f0:=f; { #### AUGMENTATION } 
   for f:=1 to f0 do 
   begin
      i1:=free[f]; low:=1; up:=1;{ ## initialize d- and pred-array } 
      for j:=1 to n do begin d[j]:=c[i1,j]-v[j]; pred[j]:=i1 end; 
      repeat 
         if up=low then  { ## find columns with new value for minimum d }
         begin 
            last:=low-1; min:=d[col[up]]; up:=up+1; 
            for k:=up to n do 
            begin 
               j:=col[k]; h:=d[j]; 
               if h<=min then 
               begin
                  if h<min then begin up:=low; min:=h end; 
                  col[k]:=col[up]; col[up]:=j; up:=up+1 
               end 
            end; 
            for h:=low to up-1 do 
            begin j:=col[h]; if y[j]=0 then goto augment end 
         end; { up=low }

         j1:=col[low]; low:=low+1; i:=y[j1];  { ## scan a row } 
         u1:=c[i,j1]-v[j1]-min; 
         for k:=up to n do 
         begin
            j:=col[k]; h:=c[i,j]-v[j]-u1; 
            if h<d[j] then 
            begin 
               d[j]:=h; pred[j]:=i; 
               if h=min then 
                   if y[j]=0 then goto augment 
                      else begin col[k]:=col[up]; col[up]:=j; up:=up+1 end 
            end 
         end 
      until false; { repeat ends with goto augment }

augment: 
     for k:=1 to last do { ## updating of column prices } 
         begin j1:=col[k]; v[j1]:=v[j1]+d[j1]-min end; 
     repeat{ ## augmentation } 
        i:=pred[j]; y[j]:=i; k:=j; j:=x[i]; x[i]:=k 
     until i=i1 
   end; { of augmentation }
 
   h:=0; { #### DETERMINE ROW PRICES AND OPTIMAL VALUE } 
   for i:=1 to n do begin j:=x[i]; u[i]:=c[i,j]-v[j]; h:=h+u[i]+v[j] end; 
   lapjv:=h 
end.
