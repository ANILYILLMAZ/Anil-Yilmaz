clc ; clear all ;
% Variables which can be changed
% Dimensions of the plate
L = 1 ;             % Length of the Plate along X-axes
B = 1 ;             % Height of the Plate along Y-axes
% Number of Elements required 

Nx = 96 ;            % Number of Elements along X-axes
Ny = Nx ;            % Number of Elements along Y-axes
l_unit=1/Nx;
%----------------------------------------
% From here dont change
nel = Nx*Ny ;        % Total Number of Elements in the Mesh
nnel = 4 ;           % Number of nodes per Element
% Number of points on the Length and Height
npx = Nx+1 ;
npy = Ny+1 ;
nnode = npx*npy ;      % Total Number of Nodes in the Mesh
% Discretizing the Length and Height of the plate
nx = linspace(0,L,npx) ;
ny = linspace(0,B,npy) ;
[xx yy] = meshgrid(nx,ny) ;
% To get the Nodal Connectivity Matrix
coordinates = [xx(:) yy(:)] ;
NodeNo = 1:nnode ;
nodes = zeros(nel,nnel) ;
% If elements along the X-axes and Y-axes are equal
% ordered counter-clockwise 
if npx==npy
    NodeNo = reshape(NodeNo,npx,npy);
    nodes(:,1) = reshape(NodeNo(1:npx-1,1:npy-1),nel,1);  %first columb of node str.
    nodes(:,2) = reshape(NodeNo(2:npx,1:npy-1),nel,1);    %second "
    nodes(:,3) = reshape(NodeNo(2:npx,2:npy),nel,1);      %thrid  "
    nodes(:,4) = reshape(NodeNo(1:npx-1,2:npy),nel,1);    %fourth  "
% If the elements along the axes are different
else%if npx>npy
    NodeNo = reshape(NodeNo,npy,npx);
    nodes(:,1) = reshape(NodeNo(1:npy-1,1:npx-1),nel,1);
    nodes(:,2) = reshape(NodeNo(2:npy,1:npx-1),nel,1);
    nodes(:,3) = reshape(NodeNo(2:npy,2:npx),nel,1);
    nodes(:,4) = reshape(NodeNo(1:npy-1,2:npx),nel,1);
end

%specified to radius of circular field % should be smaller than 0.5
radius=0.28;
%center coordinates of circular field
x0=0.5;
y0=0.5;

dx=nx(1,2)-nx(1,1);
dy=ny(1,2)-ny(1,1);

number=radius/dx;
integ=floor(number);
residual=radius-(integ*dx);

inner_radius=radius-(dx+residual*1.01);

%Make poly shape which is specified in circle in this case 

alpha_field = nsidedpoly(1000,'Center',[x0 y0],'Radius',radius);

alpha_innerfield= nsidedpoly(1000,'Center',[x0 y0],'Radius',inner_radius);


 
 
%--------------------------------------------------------------------------------%
%
% Plotting the Finite Volume Mesh
% Initialization of the required matrices
X = zeros(nnel,nel) ;
Y = zeros(nnel,nel) ;
% Extract X,Y coordinates for the (iel)-th element
  for iel = 1:nel
      X(:,iel) = coordinates(nodes(iel,:),1) ;
      Y(:,iel) = coordinates(nodes(iel,:),2) ;
  end
  



% Figure
fh = figure ;

set(fh,'name','Preprocessing for FVM','numbertitle','off','color','w') ;
patch(X,Y,'w')

title('Cartesian Grid FVM') ;
axis([0. L*1.01 0. B*1.01])
axis on ;
grid on



hold on



% To display Node Numbers % Element Numbers
pos = [70 20 60 20] ;
ShowNodes = uicontrol('style','toggle','string','nodes','value',0,....
    'position',[pos(1) pos(2) pos(3) pos(4)],'background','white','callback',...
    'SHOWNODES(ShowNodes,ShowElements,coordinates,X,Y,nnode,nel,nodes,alpha_field)');
pos = get(ShowNodes,'position') ;
pos = [2*pos(1) pos(2) pos(3) pos(4)] ;
ShowElements = uicontrol('style','toggle','string','Elements','value',0,....
    'position',[pos(1) pos(2) pos(3) pos(4)],'background','white','callback',....
    'SHOWELEMENTS(ShowElements,ShowNodes,coordinates,X,Y,nel,nodes,nnode,alpha_field)');


plot(alpha_field)
hold on
% plot(alpha_innerfield)
% hold on


el=nnode-npx-1;

inner_nodes=zeros(el,6);
%Inner node ordered with 
%[node_i,node_j,fistNeighborInXDirection_i,fistNeighborInXDirection_j,secondNeighborInYDirection_i,secondNeighborInYDirection_j]
    for nn=1:nnode
   if  ((((coordinates(nn,1)-x0)^2)+ ((coordinates(nn,2)-y0)^2))^0.5)<radius
     
    inner_nodes(nn,1)=coordinates(nn,1);
    inner_nodes(nn,2)=coordinates(nn,2);
    
    %check whether inner radius greather than inner nodes or not!
    if ((((inner_nodes(nn,1)-x0)^2)+ ((inner_nodes(nn,2)-y0)^2))^0.5)<inner_radius
    
        inner_nodes(nn,1)=0;
        inner_nodes(nn,2)=0;
    else
        
    if coordinates(nn,1)<x0 & coordinates(nn,2)<y0  %left-down portion
    
        inner_nodes(nn,3)=x0-((radius-coordinates(nn,2)+y0)*(radius+coordinates(nn,2)-y0))^(0.5);
        inner_nodes(nn,4)=coordinates(nn,2);
        inner_nodes(nn,5)=coordinates(nn,1);
        inner_nodes(nn,6)=y0-((radius-coordinates(nn,1)+y0)*(radius+coordinates(nn,1)-y0))^(0.5);
         
         
         
         
             elseif coordinates(nn,1)<x0 & coordinates(nn,2)==y0  % reduce the neighbor number on directly x and y direction 
            inner_nodes(nn,3)=x0-((radius-coordinates(nn,2)+y0)*(radius+coordinates(nn,2)-y0))^(0.5);
            inner_nodes(nn,4)=y0;  
            inner_nodes(nn,5)=coordinates(nn,1);
            inner_nodes(nn,6)=y0-((radius-coordinates(nn,1)+y0)*(radius+coordinates(nn,1)-y0))^(0.5);
            
            elseif coordinates(nn,1)==x0 & coordinates(nn,2)<y0 % reduce the neighbor number on directly x and y direction 
            inner_nodes(nn,3)=x0-((radius-coordinates(nn,2)+y0)*(radius+coordinates(nn,2)-y0))^(0.5);
            inner_nodes(nn,4)=coordinates(nn,2);    
            inner_nodes(nn,5)=coordinates(nn,1);
            inner_nodes(nn,6)=y0-((radius-coordinates(nn,1)+y0)*(radius+coordinates(nn,1)-y0))^(0.5);
            
    elseif coordinates(nn,1)<x0 & coordinates(nn,2)>y0 %left-up portion
             inner_nodes(nn,3)=x0-((radius-coordinates(nn,2)+y0)*(radius+coordinates(nn,2)-y0))^(0.5);
             inner_nodes(nn,4)=coordinates(nn,2);
             inner_nodes(nn,5)=coordinates(nn,1);
             inner_nodes(nn,6)=y0+((radius-coordinates(nn,1)+y0)*(radius+coordinates(nn,1)-y0))^(0.5);
   
    elseif coordinates(nn,1)>x0 & coordinates(nn,2)>y0 %right-up portion
             inner_nodes(nn,3)=x0+((radius-coordinates(nn,2)+y0)*(radius+coordinates(nn,2)-y0))^(0.5);
             inner_nodes(nn,4)=coordinates(nn,2);
             inner_nodes(nn,5)=coordinates(nn,1);
             inner_nodes(nn,6)=y0+((radius-coordinates(nn,1)+y0)*(radius+coordinates(nn,1)-y0))^(0.5);
             
             elseif coordinates(nn,1)==x0 & coordinates(nn,2)>y0 % reduce the neighbor number on directly x and y direction 
            inner_nodes(nn,3)=x0-((radius-coordinates(nn,2)+y0)*(radius+coordinates(nn,2)-y0))^(0.5);
            inner_nodes(nn,4)=coordinates(nn,2);    
            inner_nodes(nn,5)=coordinates(nn,1);
            inner_nodes(nn,6)=y0+((radius-coordinates(nn,1)+y0)*(radius+coordinates(nn,1)-y0))^(0.5);
             
             elseif coordinates(nn,1)>x0 & coordinates(nn,2)==y0 % reduce the neighbor number on directly x and y direction 
            inner_nodes(nn,3)=x0+((radius-coordinates(nn,2)+y0)*(radius+coordinates(nn,2)-y0))^(0.5);
            inner_nodes(nn,4)=y0;  
            inner_nodes(nn,5)=coordinates(nn,1);
            inner_nodes(nn,6)=y0+((radius-coordinates(nn,1)+y0)*(radius+coordinates(nn,1)-y0))^(0.5);
    
    elseif coordinates(nn,1)>x0 & coordinates(nn,2)<y0 %right-down portion
             inner_nodes(nn,3)=x0+((radius-coordinates(nn,2)+y0)*(radius+coordinates(nn,2)-y0))^(0.5);
             inner_nodes(nn,4)=coordinates(nn,2);
             inner_nodes(nn,5)=coordinates(nn,1);
             inner_nodes(nn,6)=y0-((radius-coordinates(nn,1)+y0)*(radius+coordinates(nn,1)-y0))^(0.5);
    
              
    else 
    end
     end
    else
    
 end 

end 
    
    
   
    
    %Re-order with Neighbour node and contact point data
    
    inner_nodes_rev=inner_nodes;
    [r,c] = size(inner_nodes);
    
    for ni=1:r
        if inner_nodes(ni,1)~=0
    %sol alt porsiyon
    %y yönü
    if inner_nodes(ni,1)<x0 & inner_nodes(ni,2)<=y0  %left-down portion
        
    if inner_nodes(ni-1,2)~=0
       
    inner_nodes_rev(ni,5)=inner_nodes(ni-1,1);
    inner_nodes_rev(ni,6)=inner_nodes(ni-1,2);
 
            
    end
    %x yönü
      if inner_nodes(ni-npx,2)~=0
       
    inner_nodes_rev(ni,3)=inner_nodes(ni-npx,1);
    inner_nodes_rev(ni,4)=inner_nodes(ni-npx,2);
 
            
    end
    
    elseif inner_nodes(ni,1)<=x0 & inner_nodes(ni,2)>y0  
    %sol üst porsiyon
    %y yönü
   
    if inner_nodes(ni+1,2)~=0
       
    inner_nodes_rev(ni,5)=inner_nodes(ni+1,1);
    inner_nodes_rev(ni,6)=inner_nodes(ni+1,2);
 
            
    end
    %x yönü
      if inner_nodes(ni-npx,2)~=0
       
    inner_nodes_rev(ni,3)=inner_nodes(ni-npx,1);
    inner_nodes_rev(ni,4)=inner_nodes(ni-npx,2);
      end
      
      
   
    
    elseif inner_nodes(ni,1)>x0 & inner_nodes(ni,2)>=y0  
    %sag üst porsiyon
    %y yönü
   
    if inner_nodes(ni+1,2)~=0
       
    inner_nodes_rev(ni,5)=inner_nodes(ni+1,1);
    inner_nodes_rev(ni,6)=inner_nodes(ni+1,2);
 
            
    end
    
    %x yönü
      if inner_nodes(ni+npx,2)~=0
       
    inner_nodes_rev(ni,3)=inner_nodes(ni+npx,1);
    inner_nodes_rev(ni,4)=inner_nodes(ni+npx,2);
 
            
      end
      
      elseif inner_nodes(ni,1)>=x0 & inner_nodes(ni,2)<y0  %left-down portion
    %sag alt porsiyon
    %y yönü
   
    if inner_nodes(ni-1,2)~=0
       
    inner_nodes_rev(ni,5)=inner_nodes(ni-1,1);
    inner_nodes_rev(ni,6)=inner_nodes(ni-1,2);
 
            
    end
    %x yönü
      if inner_nodes(ni+npx,2)~=0
       
    inner_nodes_rev(ni,3)=inner_nodes(ni+npx,1);
    inner_nodes_rev(ni,4)=inner_nodes(ni+npx,2);
 
            
      end
      
    end
        else
        end
    end
    
    
    
    %find the portional area of each active node (inner_node)
    kose=103;
    
    %dörtgenin alan?n? iki üçgenden buluyoruz ancak nokta s?ralar? ve
    %ç?karmalar?na dikkat et
    area1=abs(0.5*((inner_nodes_rev(104,1)*(inner_nodes_rev(103,2)-inner_nodes_rev(103,4))+inner_nodes_rev(103,1)*(inner_nodes_rev(103,4)-inner_nodes_rev(104,2))+inner_nodes_rev(103,3)*(inner_nodes_rev(104,2)-inner_nodes_rev(103,2)))));
    area2=abs(0.5*((inner_nodes_rev(104,1)*(inner_nodes_rev(103,4)-inner_nodes_rev(104,4))+inner_nodes_rev(103,3)*(inner_nodes_rev(104,4)-inner_nodes_rev(104,2))+inner_nodes_rev(104,3)*(inner_nodes_rev(104,2)-inner_nodes_rev(103,4)))));
    
    area=area1+area2;
    general_area=dx*dy;
    portional_area=area/general_area;
    
    
    
    %must order with clock-wise direction
    
    for ar=npy+1:r
    if inner_nodes_rev(ar,1)<x0 & inner_nodes_rev(ar,2)<y0  %left-down portion
        
    if ((inner_nodes_rev(ar,3)==inner_nodes_rev(ar-npx,1))&(inner_nodes_rev(ar,4)==inner_nodes_rev(ar-npx,2))) & ((inner_nodes_rev(ar,5)==inner_nodes_rev(ar-1,1))&(inner_nodes_rev(ar,6)==inner_nodes_rev(ar-1,2)))
    % 5 side polygon area calculation
  
    P = [inner_nodes_rev(ar,1) inner_nodes_rev(ar,2); inner_nodes_rev(ar,5) inner_nodes_rev(ar,6);inner_nodes_rev(ar-1,3) inner_nodes_rev(ar-1,4);inner_nodes_rev(ar-npx,5) inner_nodes_rev(ar-npx,6);inner_nodes_rev(ar,3) inner_nodes_rev(ar,4);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    
    A(ar,1)=polyin.area/general_area;
    
    elseif ((inner_nodes_rev(ar,3)==inner_nodes_rev(ar-npx,1))&(inner_nodes_rev(ar,4)==inner_nodes_rev(ar-npx,2))) & ((inner_nodes_rev(ar,5)~=inner_nodes_rev(ar-1,1))&(inner_nodes_rev(ar,6)~=inner_nodes_rev(ar-1,2))) | ((inner_nodes_rev(ar,3)~=inner_nodes_rev(ar-npx,1))&(inner_nodes_rev(ar,4)~=inner_nodes_rev(ar-npx,2))) & ((inner_nodes_rev(ar,5)==inner_nodes_rev(ar-1,1))&(inner_nodes_rev(ar,6)==inner_nodes_rev(ar-1,2)))
    % 4 side polygon area calculation
    %sol alt kose ust porsiyon
    if inner_nodes_rev(ar,2) > inner_nodes_rev(ar,1)
    P = [ inner_nodes_rev(ar,3) inner_nodes_rev(ar,4); inner_nodes_rev(ar,1) inner_nodes_rev(ar,2);inner_nodes_rev(ar,5) inner_nodes_rev(ar,6); inner_nodes_rev(ar-1,3) inner_nodes_rev(ar-1,4);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    A(ar,1)=polyin.area/general_area;
    
    else %sol alt kose alt porsiyon
    P = [ inner_nodes_rev(ar,3) inner_nodes_rev(ar,4); inner_nodes_rev(ar,1) inner_nodes_rev(ar,2);inner_nodes_rev(ar,5) inner_nodes_rev(ar,6); inner_nodes_rev(ar-npx,5) inner_nodes_rev(ar-npx,6);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    A(ar,1)=polyin.area/general_area;
    
    
    end
    elseif ((inner_nodes_rev(ar,3)~=inner_nodes_rev(ar-npx,1))&(inner_nodes_rev(ar,4)~=inner_nodes_rev(ar-npx,2))) & ((inner_nodes_rev(ar,5)~=inner_nodes_rev(ar-1,1))&(inner_nodes_rev(ar,6)~=inner_nodes_rev(ar-1,2)))
    % 3 side triangle area calculation
         P = [ inner_nodes_rev(ar,3) inner_nodes_rev(ar,4); inner_nodes_rev(ar,1) inner_nodes_rev(ar,2);inner_nodes_rev(ar,5) inner_nodes_rev(ar,6);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    A(ar,1)=polyin.area/general_area;
        
    else
        cevap="yanlis"
        
       
        
        
    end
    
    
    elseif inner_nodes_rev(ar,1)<x0 & inner_nodes_rev(ar,2)>y0  %left-up portion
         
        
    if ((inner_nodes_rev(ar,3)==inner_nodes_rev(ar-npx,1))&(inner_nodes_rev(ar,4)==inner_nodes_rev(ar-npx,2))) & ((inner_nodes_rev(ar,5)==inner_nodes_rev(ar+1,1))&(inner_nodes_rev(ar,6)==inner_nodes_rev(ar+1,2)))
    % 5 side polygon area calculation
    P = [inner_nodes_rev(ar,1) inner_nodes_rev(ar,2); inner_nodes_rev(ar,3) inner_nodes_rev(ar,4);inner_nodes_rev(ar-npx,5) inner_nodes_rev(ar-npx,6);inner_nodes_rev(ar+1,3) inner_nodes_rev(ar+1,4);inner_nodes_rev(ar,5) inner_nodes_rev(ar,6);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    
    A(ar,1)=polyin.area/general_area;
    
    elseif ((inner_nodes_rev(ar,3)==inner_nodes_rev(ar-npx,1))&(inner_nodes_rev(ar,4)==inner_nodes_rev(ar-npx,2))) & ((inner_nodes_rev(ar,5)~=inner_nodes_rev(ar+1,1))&(inner_nodes_rev(ar,6)~=inner_nodes_rev(ar+1,2))) | ((inner_nodes_rev(ar,3)~=inner_nodes_rev(ar-npx,1))&(inner_nodes_rev(ar,4)~=inner_nodes_rev(ar-npx,2))) & ((inner_nodes_rev(ar,5)==inner_nodes_rev(ar+1,1))&(inner_nodes_rev(ar,6)==inner_nodes_rev(ar+1,2)))
    % 4 side polygon area calculation
    
    if inner_nodes_rev(ar,2)>(2*y0)-inner_nodes_rev(ar,1)
    %sol-ust ust porsiyon
    P = [ inner_nodes_rev(ar,3) inner_nodes_rev(ar,4); inner_nodes_rev(ar-npx,5) inner_nodes_rev(ar-npx,6);inner_nodes_rev(ar,5) inner_nodes_rev(ar,6); inner_nodes_rev(ar,1) inner_nodes_rev(ar,2);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    A(ar,1)=polyin.area/general_area;
    
    else
        
    %sol-ust alt porsiyon
    P = [ inner_nodes_rev(ar,3) inner_nodes_rev(ar,4); inner_nodes_rev(ar+1,3) inner_nodes_rev(ar+1,4);inner_nodes_rev(ar,5) inner_nodes_rev(ar,6); inner_nodes_rev(ar,1) inner_nodes_rev(ar,2);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    A(ar,1)=polyin.area/general_area ;
        
        
        
    end
    elseif ((inner_nodes_rev(ar,3)~=inner_nodes_rev(ar-npx,1))&(inner_nodes_rev(ar,4)~=inner_nodes_rev(ar-npx,2))) & ((inner_nodes_rev(ar,5)~=inner_nodes_rev(ar+1,1))&(inner_nodes_rev(ar,6)~=inner_nodes_rev(ar+1,2)))
    % 3 side triangle area calculation
         P = [ inner_nodes_rev(ar,3) inner_nodes_rev(ar,4); inner_nodes_rev(ar,1) inner_nodes_rev(ar,2);inner_nodes_rev(ar,5) inner_nodes_rev(ar,6);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    A(ar,1)=polyin.area/general_area;
        
    else
        cevap="yanlis"
        
       
        
        
    end
        
        
        elseif inner_nodes_rev(ar,1)>x0 & inner_nodes_rev(ar,2)>y0  %right-up portion
            
              
    if ((inner_nodes_rev(ar,3)==inner_nodes_rev(ar+npx,1))&(inner_nodes_rev(ar,4)==inner_nodes_rev(ar+npx,2))) & ((inner_nodes_rev(ar,5)==inner_nodes_rev(ar+1,1))&(inner_nodes_rev(ar,6)==inner_nodes_rev(ar+1,2)))
    % 5 side polygon area calculation
    P = [inner_nodes_rev(ar,1) inner_nodes_rev(ar,2); inner_nodes_rev(ar,3) inner_nodes_rev(ar,4);inner_nodes_rev(ar+npx,5) inner_nodes_rev(ar+npx,6);inner_nodes_rev(ar+1,3) inner_nodes_rev(ar+1,4);inner_nodes_rev(ar,5) inner_nodes_rev(ar,6);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    
    A(ar,1)=polyin.area/general_area;
    
    elseif ((inner_nodes_rev(ar,3)==inner_nodes_rev(ar+npx,1))&(inner_nodes_rev(ar,4)==inner_nodes_rev(ar+npx,2))) & ((inner_nodes_rev(ar,5)~=inner_nodes_rev(ar+1,1))&(inner_nodes_rev(ar,6)~=inner_nodes_rev(ar+1,2))) | ((inner_nodes_rev(ar,3)~=inner_nodes_rev(ar+npx,1))&(inner_nodes_rev(ar,4)~=inner_nodes_rev(ar+npx,2))) & ((inner_nodes_rev(ar,5)==inner_nodes_rev(ar+1,1))&(inner_nodes_rev(ar,6)==inner_nodes_rev(ar+1,2)))
    % 4 side polygon area calculation
    
    if inner_nodes_rev(ar,2)>inner_nodes_rev(ar,1)
    %sol-ust ust porsiyon
    P = [ inner_nodes_rev(ar,3) inner_nodes_rev(ar,4); inner_nodes_rev(ar+npx,5) inner_nodes_rev(ar+npx,6);inner_nodes_rev(ar,5) inner_nodes_rev(ar,6); inner_nodes_rev(ar,1) inner_nodes_rev(ar,2);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    A(ar,1)=polyin.area/general_area;
    
    else
        
    %sol-ust alt porsiyon
    P = [ inner_nodes_rev(ar,3) inner_nodes_rev(ar,4); inner_nodes_rev(ar+1,3) inner_nodes_rev(ar+1,4);inner_nodes_rev(ar,5) inner_nodes_rev(ar,6); inner_nodes_rev(ar,1) inner_nodes_rev(ar,2);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    A(ar,1)=polyin.area/general_area ;
        
        
        
    end
    elseif ((inner_nodes_rev(ar,3)~=inner_nodes_rev(ar+npx,1))&(inner_nodes_rev(ar,4)~=inner_nodes_rev(ar+npx,2))) & ((inner_nodes_rev(ar,5)~=inner_nodes_rev(ar+1,1))&(inner_nodes_rev(ar,6)~=inner_nodes_rev(ar+1,2)))
    % 3 side triangle area calculation
         P = [ inner_nodes_rev(ar,3) inner_nodes_rev(ar,4); inner_nodes_rev(ar,1) inner_nodes_rev(ar,2);inner_nodes_rev(ar,5) inner_nodes_rev(ar,6);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    A(ar,1)=polyin.area/general_area;
        
    else
        cevap="yanlis"
        
       
        
        
    end  
    elseif inner_nodes_rev(ar,1)>x0 & inner_nodes_rev(ar,2)<y0  %left-down portion
        
    if ((inner_nodes_rev(ar,3)==inner_nodes_rev(ar+npx,1))&(inner_nodes_rev(ar,4)==inner_nodes_rev(ar+npx,2))) & ((inner_nodes_rev(ar,5)==inner_nodes_rev(ar-1,1))&(inner_nodes_rev(ar,6)==inner_nodes_rev(ar-1,2)))
    % 5 side polygon area calculation
  
    P = [inner_nodes_rev(ar,1) inner_nodes_rev(ar,2); inner_nodes_rev(ar,5) inner_nodes_rev(ar,6);inner_nodes_rev(ar-1,3) inner_nodes_rev(ar-1,4);inner_nodes_rev(ar+npx,5) inner_nodes_rev(ar+npx,6);inner_nodes_rev(ar,3) inner_nodes_rev(ar,4);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    
    A(ar,1)=polyin.area/general_area;
    
    elseif ((inner_nodes_rev(ar,3)==inner_nodes_rev(ar+npx,1))&(inner_nodes_rev(ar,4)==inner_nodes_rev(ar+npx,2))) & ((inner_nodes_rev(ar,5)~=inner_nodes_rev(ar-1,1))&(inner_nodes_rev(ar,6)~=inner_nodes_rev(ar-1,2))) | ((inner_nodes_rev(ar,3)~=inner_nodes_rev(ar+npx,1))&(inner_nodes_rev(ar,4)~=inner_nodes_rev(ar+npx,2))) & ((inner_nodes_rev(ar,5)==inner_nodes_rev(ar-1,1))&(inner_nodes_rev(ar,6)==inner_nodes_rev(ar-1,2)))
    % 4 side polygon area calculation
    %sol alt kose ust porsiyon
    if inner_nodes_rev(ar,2) >(2*y0)- inner_nodes_rev(ar,1)
    P = [ inner_nodes_rev(ar,3) inner_nodes_rev(ar,4); inner_nodes_rev(ar,1) inner_nodes_rev(ar,2);inner_nodes_rev(ar,5) inner_nodes_rev(ar,6); inner_nodes_rev(ar-1,3) inner_nodes_rev(ar-1,4);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    A(ar,1)=polyin.area/general_area;
    
    else %sol alt kose alt porsiyon
    P = [ inner_nodes_rev(ar,3) inner_nodes_rev(ar,4); inner_nodes_rev(ar,1) inner_nodes_rev(ar,2);inner_nodes_rev(ar,5) inner_nodes_rev(ar,6); inner_nodes_rev(ar+npx,5) inner_nodes_rev(ar+npx,6);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    A(ar,1)=polyin.area/general_area;
    
    
    end
    elseif ((inner_nodes_rev(ar,3)~=inner_nodes_rev(ar+npx,1))&(inner_nodes_rev(ar,4)~=inner_nodes_rev(ar+npx,2))) & ((inner_nodes_rev(ar,5)~=inner_nodes_rev(ar-1,1))&(inner_nodes_rev(ar,6)~=inner_nodes_rev(ar-1,2)))
    % 3 side triangle area calculation
         P = [ inner_nodes_rev(ar,3) inner_nodes_rev(ar,4); inner_nodes_rev(ar,1) inner_nodes_rev(ar,2);inner_nodes_rev(ar,5) inner_nodes_rev(ar,6);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    A(ar,1)=polyin.area/general_area;
        
    else
        cevap="yanlis"
        
       
        
        
    end
     
    elseif inner_nodes_rev(ar,1)<x0 & inner_nodes_rev(ar,2)==y0 
        
    if ((inner_nodes_rev(ar,3)==inner_nodes_rev(ar-npx,1))&(inner_nodes_rev(ar,4)==inner_nodes_rev(ar-npx,2))) & ((inner_nodes_rev(ar,5)==inner_nodes_rev(ar-1,1))&(inner_nodes_rev(ar,6)==inner_nodes_rev(ar-1,2)))
    % 5 side polygon area calculation
  
    P = [inner_nodes_rev(ar,1) inner_nodes_rev(ar,2); inner_nodes_rev(ar-1,1) inner_nodes_rev(ar-1,2);inner_nodes_rev(ar-1,3) inner_nodes_rev(ar-1,4);inner_nodes_rev(ar-npx,5) inner_nodes_rev(ar-npx,6);inner_nodes_rev(ar,3) inner_nodes_rev(ar,4);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    
       
    
    A(ar,1)=polyin.area/general_area;
    
    girdi_alan=A(ar,1)
    
    elseif ((inner_nodes_rev(ar,3)==inner_nodes_rev(ar-npx,1))&(inner_nodes_rev(ar,4)==inner_nodes_rev(ar-npx,2))) & ((inner_nodes_rev(ar,5)~=inner_nodes_rev(ar-1,1))&(inner_nodes_rev(ar,6)~=inner_nodes_rev(ar-1,2))) | ((inner_nodes_rev(ar,3)~=inner_nodes_rev(ar-npx,1))&(inner_nodes_rev(ar,4)~=inner_nodes_rev(ar-npx,2))) & ((inner_nodes_rev(ar,5)==inner_nodes_rev(ar-1,1))&(inner_nodes_rev(ar,6)==inner_nodes_rev(ar-1,2)))
    % 4 side polygon area calculation
    %sol alt kose ust porsiyon
    if inner_nodes_rev(ar,2) > inner_nodes_rev(ar,1)
    P = [ inner_nodes_rev(ar,3) inner_nodes_rev(ar,4); inner_nodes_rev(ar,1) inner_nodes_rev(ar,2);inner_nodes_rev(ar,5) inner_nodes_rev(ar,6); inner_nodes_rev(ar-1,3) inner_nodes_rev(ar-1,4);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    A(ar,1)=polyin.area/general_area;
    
    else %sol alt kose alt porsiyon
    P = [ inner_nodes_rev(ar,3) inner_nodes_rev(ar,4); inner_nodes_rev(ar,1) inner_nodes_rev(ar,2);inner_nodes_rev(ar,5) inner_nodes_rev(ar,6); inner_nodes_rev(ar-npx,5) inner_nodes_rev(ar-npx,6);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    A(ar,1)=polyin.area/general_area;
    
    
    end
    elseif ((inner_nodes_rev(ar,3)~=inner_nodes_rev(ar-npx,1))&(inner_nodes_rev(ar,4)~=inner_nodes_rev(ar-npx,2))) & ((inner_nodes_rev(ar,5)~=inner_nodes_rev(ar-1,1))&(inner_nodes_rev(ar,6)~=inner_nodes_rev(ar-1,2)))
    % 3 side triangle area calculation
         P = [ inner_nodes_rev(ar,3) inner_nodes_rev(ar,4); inner_nodes_rev(ar,1) inner_nodes_rev(ar,2);inner_nodes_rev(ar,5) inner_nodes_rev(ar,6);];

    polyin= polyshape(P);
    
    plot(polyin)
    axis equal
    A(ar,1)=polyin.area/general_area;
        
    else
        cevap="yanlis"
        
       
        
        
    end
    
       
    end
    end
    
    
%Distrubition of volume fraction (alpha) on each nodes and elements

alpha=zeros(nel,1);
alpha(:,1)=0;





%Make distrubition matrices that includes volume fraction field
dist=[nodes(:,1:2:3:4)]; 
     

    
   % collect area information and element number 
   for sycc=1:ni
       
   if inner_nodes_rev(sycc,1)~=0 & inner_nodes_rev(sycc,1)<=x0 & inner_nodes_rev(sycc,2)<=y0
       
   for sycel=1:nel
   if dist(sycel,3)==sycc 
       
       alpha(sycel,1)=A(sycc,1);
       
       
   end
   
   end
   
   elseif inner_nodes_rev(sycc,1)~=0 & inner_nodes_rev(sycc,1)<=x0 & inner_nodes_rev(sycc,2)>=y0
       
   for sycel=1:nel
   if dist(sycel,4)==sycc 
       
       alpha(sycel,1)=A(sycc,1);
       
       
   end
   
   end
   elseif inner_nodes_rev(sycc,1)~=0 & inner_nodes_rev(sycc,1)>=x0 & inner_nodes_rev(sycc,2)<=y0
       
   for sycel=1:nel
   if dist(sycel,2)==sycc 
       
       alpha(sycel,1)=A(sycc,1);
       
       
   end
   
   end
       
        elseif inner_nodes_rev(sycc,1)~=0 & inner_nodes_rev(sycc,1)>=x0 & inner_nodes_rev(sycc,2)>=y0
       
   for sycel=1:nel
   if dist(sycel,1)==sycc 
       
       alpha(sycel,1)=A(sycc,1);
       
       
   end
   
   end
   end
   end
   
   gg=1;
   
   %re-order with axial data x and y
   
   for sycc=1:ni
       
   if  inner_nodes_rev(sycc,2)==y0 & inner_nodes_rev(sycc,1)<x0
       
   for sycel=1:nel
   if dist(sycel,3)==sycc 
       
       alpha(sycel,1)=A(sycc,1);
       yakala(gg)=sycel;  %keep the same area value for put to other 7 axial area
       gg=gg+1;
   elseif dist(sycel,4)==sycc 
       
       alpha(sycel,1)=alpha(yakala(1),1);
       
   end
   
   end
    
   elseif  inner_nodes_rev(sycc,2)==y0 & inner_nodes_rev(sycc,1)>x0
       
   for sycel=1:nel
   if dist(sycel,2)==sycc 
       
       alpha(sycel,1)=alpha(yakala(1),1);
       
   elseif dist(sycel,1)==sycc 
       
       alpha(sycel,1)=alpha(yakala(1),1);
   end
   
   end
   elseif  inner_nodes_rev(sycc,1)==x0 & inner_nodes_rev(sycc,2)>y0
       
   for sycel=1:nel
   if dist(sycel,4)==sycc 
       
       alpha(sycel,1)=alpha(yakala(1),1);
       
   elseif dist(sycel,1)==sycc 
       
       alpha(sycel,1)=alpha(yakala(1),1);
       
       
   end
   
   end
   elseif  inner_nodes_rev(sycc,1)==x0 & inner_nodes_rev(sycc,2)<y0
       
   for sycel=1:nel
   if dist(sycel,3)==sycc 
       
       alpha(sycel,1)=alpha(yakala(1),1);
       
   elseif dist(sycel,2)==sycc 
       
       alpha(sycel,1)=alpha(yakala(1),1);
       
       
   end
   
   end
   
   end
   end
   
%    filled into circle element those nodes are fully in!!  "alpha=1"
   
   for i_el=1:nel
       
       %1,5,9,13
       
      sec(1,1)=dist(i_el,1);
      sec(2,1)=dist(i_el,2);
      sec(3,1)=dist(i_el,3);
      sec(4,1)=dist(i_el,4);
      
      if ((((coordinates(sec(1,1),1)-x0)^2)+ ((coordinates(sec(1,1),2)-y0)^2))^0.5)<=radius &((((coordinates(sec(2,1),1)-x0)^2)+ ((coordinates(sec(2,1),2)-y0)^2))^0.5)<=radius &((((coordinates(sec(3,1),1)-x0)^2)+ ((coordinates(sec(3,1),2)-y0)^2))^0.5)<=radius &((((coordinates(sec(4,1),1)-x0)^2)+ ((coordinates(sec(4,1),2)-y0)^2))^0.5)<=radius
        
       alpha(i_el,1)=1;
      end
       
       
       
       
       
   end
   
   
   
   
   
 

snormal=zeros(nel,8);
snormal(:,1)=-1;
snormal(:,2)=0;
snormal(:,3)=0;
snormal(:,4)=1;
snormal(:,5)=1;
snormal(:,6)=0;
snormal(:,7)=0;
snormal(:,8)=-1;

%Surface normal's and between two nodes surfaces
%Ordered: ( (nodes ,i,j,nodes) .... X4 ,first nodes, volume fraction) /// direction : counter-clockwise  

dist_snormal_alpha=[dist(:,1) snormal(:,1:2) dist(:,2) snormal(:,3:4) dist(:,3) snormal(:,5:6) dist(:,4) snormal(:,7:8) dist(:,1) alpha  ];

s_value=zeros(nel,4);
lim_u=Nx*Ny-(Nx+1);
lim_l= Nx+2;

%surface value calculation
 for el = lim_l:lim_u
  
   if el==Nx*2 || el==1+Nx*2 || el==Nx*3 || el==1+Nx*3 || el==Nx*4 || el==1+Nx*4 || el==Nx*5 || el==1+Nx*5 || el==Nx*6 || el==1+Nx*6 || el==Nx*7 || el==1+Nx*7  || el==Nx*8 || el==1+Nx*8 || el==Nx*9 || el==1+Nx*9 || el==Nx*10 || el==1+Nx*10  || el==Nx*11 || el==1+Nx*11 || el==Nx*12 || el==1+Nx*12 || el==Nx*13 || el==1+Nx*13 || el==Nx*14 || el==1+Nx*14 || el==Nx*15 || el==1+Nx*15  || el==Nx*16 || el==1+Nx*16  
   
   else
   mid=el;
   midl=el-Nx;
   midu=el+1;
   midr=el+Nx;
   midd=el-1;

s_value(mid,1)=(dist_snormal_alpha(mid,14)+dist_snormal_alpha(midl,14))/2;
s_value(mid,2)=(dist_snormal_alpha(mid,14)+dist_snormal_alpha(midu,14))/2;
s_value(mid,3)=(dist_snormal_alpha(mid,14)+dist_snormal_alpha(midr,14))/2;
s_value(mid,4)=(dist_snormal_alpha(mid,14)+dist_snormal_alpha(midd,14))/2;   
       
   end
       
  

 end

 %Distributing all variable //node number // surface normal ..
 dist_snormal_alpha=[dist(:,1) snormal(:,1:2) s_value(:,1) dist(:,2) snormal(:,3:4) s_value(:,2) dist(:,3) snormal(:,5:6) s_value(:,3) dist(:,4) snormal(:,7:8) s_value(:,4 ) dist(:,1) alpha  ];
   
   
   
   
   

 
 %Gradient calculation 

 %attention! there are no area and volume value, should add these variable
 
 for n = 1:nel
      
  gradient_of_element(n,1)=  (dist_snormal_alpha(n,4)*dist_snormal_alpha(n,2)) + (dist_snormal_alpha(n,12)*dist_snormal_alpha(n,10)) ;
 
  gradient_of_element(n,2)=  (dist_snormal_alpha(n,8)*dist_snormal_alpha(n,7)) + (dist_snormal_alpha(n,16)*dist_snormal_alpha(n,15));

  gradient_of_element(n,3)= (gradient_of_element(n,1)^2 +gradient_of_element(n,2)^2)^0.5 ;

 end
 
 %find the local centeral coordinates and put grad matrix 4 and 5 columns
 for n_el=1:nel
       
     sec_g(1,1)=dist(n_el,1);
     sec_g(2,1)=dist(n_el,2);
     sec_g(3,1)=dist(n_el,3);
     sec_g(4,1)=dist(n_el,4);
       
     gradient_of_element(n_el,4)=(coordinates(sec_g(1,1),1)+coordinates(sec_g(2,1),1)+coordinates(sec_g(3,1),1)+coordinates(sec_g(4,1),1))/4;
     gradient_of_element(n_el,5)=(coordinates(sec_g(1,1),2)+coordinates(sec_g(2,1),2)+coordinates(sec_g(3,1),2)+coordinates(sec_g(4,1),2))/4;
      
      
       
       
       
 end

 
%find max value of gradient vector
[row, col] = find(ismember(round(gradient_of_element(:,3),3), max(round(gradient_of_element(:,3),3))));
 
 
 
     
figure
plot(alpha_field)
hold on
quiverwcolorbar(gradient_of_element(:,4),gradient_of_element(:,5),gradient_of_element(:,1),gradient_of_element(:,2),0.15,jet) % scale factor must be 0-0.2 range
% colormap jet
% colorbar
title (['Gradient Distribution - Unit Lenght Of Elements:',num2str(l_unit)]);

figure
plot(alpha_field)
hold on
quiver(gradient_of_element(row(:),4),gradient_of_element(row(:),5),gradient_of_element(row(:),1),gradient_of_element(row(:),2),0.5,'r')
title ('Maximum Gradient Vector Position')

 

for n=1:iel
    
      if gradient_of_element(n,4)>=x0 && gradient_of_element(n,5)>=y0
    thetas(n,1)=atand(abs(gradient_of_element(n,5)-y0)/(abs(gradient_of_element(n,4)-x0)));
    
    elseif gradient_of_element(n,4)<x0 && gradient_of_element(n,5)>y0
        
        thetas(n,1)=180-atand(abs(gradient_of_element(n,5)-y0)/(abs(gradient_of_element(n,4)-x0)));
        
    elseif gradient_of_element(n,4)<x0 && gradient_of_element(n,5)<y0
        
                thetas(n,1)=180+atand(abs(gradient_of_element(n,5)-y0)/(abs(gradient_of_element(n,4)-x0)));
                
                
    elseif gradient_of_element(n,4)>x0 && gradient_of_element(n,5)<y0
        
            thetas(n,1)=360-atand(abs(gradient_of_element(n,5)-y0)/(abs(gradient_of_element(n,4)-x0)));

    end

   
%     thetas(n,1)=atand((abs(gradient_of_element(n,5))-abs(y0))/(abs(gradient_of_element(n,4))-abs(x0)));
    r_dist(n,1)=((gradient_of_element(n,5)-y0)^2 + (gradient_of_element(n,4)-x0)^2)^0.5;
    
end


gradient_of_element(:,6)=thetas(:,1);  % theta degree of each element

gradient_of_element(:,7)=r_dist(:,1); % centeral distance of each element

figure
scatter(gradient_of_element(:,6),gradient_of_element(:,3),'filled');
xlabel('Theta Degree')
xticks(0:60:360)
ylabel('Magnitude Of Gradient Vector')
yticks(0:0.1:1)
gradient_normalized(:,1)=gradient_of_element(:,3)/gradient_of_element(row(1,1),3);
gradient_normalized(:,2)= gradient_of_element(:,6);
figure
scatter(gradient_normalized(:,2),gradient_normalized(:,1),'filled');
xlabel('Theta Degree')
xticks(0:45:360)
ylabel('Magnitude Of Gradient Vector_Normalized')
yticks(0:0.1:1)
x0=700;
y0=300;
width=625;
height=500
set(gcf,'position',[x0,y0,width,height])