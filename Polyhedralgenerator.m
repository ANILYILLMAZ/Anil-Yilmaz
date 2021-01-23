clc ; clear all ;
% SECTION TITLE
% DESCRIPTIVE TEXT
% Variables which can be changed
% Number of Elements required 

Nx =48;            % Number of Elements along X-axes *** should take only even numbers 
Ny = Nx ;            % Number of Elements along Y-axes

% Dimensions of the plate
%L = 6.92820323 ;             % Length of the Plate along X-axes
L = 1 ; 
% From here dont change
dx=L/Nx;
dy=(dx*(3^(0.5))/2);
B =((Nx-1)*dy)+(dx*2/3^(0.5));             % Height of the Plate along Y-axes
l_unit=L/(Nx*(3^0.5));   %unit side lenght of each element

general_area=((6*l_unit^2)*3^0.5)/4;
%----------------------------------------

%Distrubition of volume fraction (alpha) on each nodes and elements

alpha=zeros((Nx^2)-(Nx/2),1);
%specified to radius of circular field % should be smaller than 0.5
radius=0.25;

x0=0.5;
y0=0.5132;
resolution=1000;

number=radius/dx;
integ=floor(number);
residual=radius-(integ*dx);

inner_radius=radius-(2*l_unit);
%Make poly shape which is specified in circle in this case 

alpha_field = nsidedpoly(resolution,'Center',[x0 y0],'Radius',radius);

alpha_inner_field = nsidedpoly(resolution,'Center',[x0 y0],'Radius',inner_radius);

%----------------------------------------
% From here dont change
nel = Nx*Ny ;        % Total Number of Elements in the Mesh
nnel = 6 ;           % Number of nodes per Element
% Number of points on the Length and Height
npx = (2*Nx)+1 ;
npy = (3*Ny)+2 ;
nnode = npx*npy ;      % Total Number of Nodes in the Mesh


% Discretizing the Length and Height of the plate
nx = linspace(0,L,npx) ;
ny = linspace(0,B,npy) ;
[xx yy] = meshgrid(nx,ny) ;

hh=1;
for syc =1:Nx/2
    a(syc,1)=hh;
    hh=hh+6;
end
%Eliminate un-useful nodes according to even or not numbers
for j=1:npx
       if (mod(j,2)==0)
           xx(2,j) =0;
               for syc2=1:Nx/2
            for i=[2,3,5,7]
           xx(a(syc2,1)+i,j) =0;   
            end
               end
       else
            xx(1,j)=0;
             for syc2=1:Nx/2
            for i=[2,4,5,6]
            xx(a(syc2,1)+i,j)=0;
            end
             end
        end
end

for j=1:npx
       if (mod(j,2)==0)
           yy(2,j) =0;
           for syc2=1:Nx/2 
           for i=[2,3,5,7]
           yy(a(syc2,1)+i,j) =0;   
           end
           end
       else
           yy(1,j) =0;
           for syc2=1:Nx/2
           for i=[2,4,5,6]
            yy(a(syc2,1)+i,j)=0;
           end
           end
        end
 end
 
%{
figure(1)
hold on
for j=1:npx
    plot(xx(:,j),yy(:,j),'k*')
end
 
for i=1:npy
    plot(xx(i,:),yy(i,:),'k*')
end
 
xlabel('x','FontSize',16)
ylabel('y','FontSize',16)
set(gca,'FontSize',16)
hold off

%}
% To get the Nodal Connectivity Matrix
coordinates = [xx(:) yy(:)] ;

%Eliminate to all zero rows in coordinate matrix
revcoordinates = coordinates;
revcoordinates( all(~revcoordinates,2), : ) = [];  % removes all rows with all zero

[r,c] = size(revcoordinates);

NodeNo = 1:r ;

%find the lateral node number in one x section 
n=((2*Nx)+1);
%find the node number in one y section 
s = r/n;
nodes = zeros((Nx^2)-(Nx/2),nnel) ;
% If elements along the X-axes and Y-axes are equal ordered
% counter-clockwise

    NodeNo = reshape(NodeNo,s,n);
   
    counter=1;
    for nu=1:Nx
        
             num=2;
        if (mod(nu,2)==0)
            for t=1:Nx-1
    nodes(counter,1) = NodeNo(nu,num);
    nodes(counter,2) = NodeNo(nu+1,num);
    nodes(counter,3) = NodeNo(nu+1,num+1);
    nodes(counter,4) = NodeNo(nu+1,num+2);
    nodes(counter,5) = NodeNo(nu,num+2);
    nodes(counter,6) = NodeNo(nu,num+1);
    counter=counter+1;
    num=num+2;
            end
        else 
             num=1;
        for t=1:Nx
            
    nodes(counter,1) = NodeNo(nu,num);
    nodes(counter,2) = NodeNo(nu+1,num);
    nodes(counter,3) = NodeNo(nu+1,num+1);
    nodes(counter,4) = NodeNo(nu+1,num+2);
    nodes(counter,5) = NodeNo(nu,num+2);
    nodes(counter,6) = NodeNo(nu,num+1);
    counter=counter+1;
    num=num+2;
        end
      end
    end





%Make distrubition matrices that includes volume fraction field 
%oriented with node number [1,2,3,4,5,6,element_alpha_fraction]
dist=[nodes(:,1:2:3:4:5:6) alpha]; 

snormal=zeros((Nx^2)-(Nx/2),8);
snormal(:,1)=-1;
snormal(:,2)=0;
snormal(:,3)=-1/2;
snormal(:,4)=(3^(0.5))/2;
snormal(:,5)=1/2;
snormal(:,6)=(3^(0.5))/2;
snormal(:,7)=1;
snormal(:,8)=0;
snormal(:,9)=1/2;
snormal(:,10)=-(3^(0.5))/2;
snormal(:,11)=-1/2;
snormal(:,12)=-(3^(0.5))/2;





%Surface normal's and between two nodes surfaces
%Ordered: ( (nodes ,i,j,nodes) .... X6 ,first nodes, volume fraction) /// direction : counter-clockwise  

dist_snormal_alpha=[dist(:,1) snormal(:,1:2) dist(:,2) snormal(:,3:4) dist(:,3) snormal(:,5:6) dist(:,4) snormal(:,7:8) dist(:,5) snormal(:,9:10) dist(:,6) snormal(:,11:12) dist(:,1) alpha  ];

%!!!!!!!!!!!!!!!! 12/04/2020 en son burada kaldim!!!!!!!!!!!!!!!!!!!!111


%This calculation for extract to exist element (exist all neighbor, especiaaly centeral area of distribution) on gradient computation  
s_value=zeros((Nx^2)-(Nx/2),6);
lim_u=(Nx^2)-(Nx/2)-Nx;
lim_l= Nx+2;

n_active_el=(Nx^2)-(Nx/2);

n_counter=(Nx-3);


%nbos_el means that extract to omitted cell number into a grid
nbos_el(1,1)= n_active_el/(Nx/2);
nbos_el(1,2)= nbos_el(1,1)+1;

for sayac=2:1:n_counter
    
 if (mod(sayac,2)==0)
nbos_el(sayac,1)= nbos_el(sayac-1,1)+Nx;
nbos_el(sayac,2)= nbos_el(sayac-1,2)+Nx;
 
 else
 nbos_el(sayac,1)= nbos_el(sayac-1,1)+Nx-1;
 nbos_el(sayac,2)= nbos_el(sayac-1,2)+Nx-1;
 end
end

nbos_el=reshape(nbos_el',[],1);




 
 
 
 
 
%--------------------------------------------------------------------------------%
%
% Plotting the Finite Volume Mesh
% Initialization of the required matrices
X = zeros(nnel,(Nx^2)-(Nx/2)) ;
Y = zeros(nnel,(Nx^2)-(Nx/2)) ;
% Extract X,Y coordinates for the (iel)-th element
  for iel = 1:(Nx^2)-(Nx/2)
      X(:,iel) = revcoordinates(nodes(iel,:),1) ;
      Y(:,iel) = revcoordinates(nodes(iel,:),2) ;
  end
  


% Figure
fh = figure ;

set(fh,'name','Preprocessing for FVM','numbertitle','off','color','w') ;
patch(X,Y,'w')
title('Polyhedral Grid FVM') ;
axis([0. L*1.01 0. B*1.01])
axis on ;
grid on

if L==B
    axis equal ;
end
hold on





% To display Node Numbers % Element Numbers
pos = [70 20 60 20] ;
ShowNodes = uicontrol('style','toggle','string','nodes','value',0,....
    'position',[pos(1) pos(2) pos(3) pos(4)],'background','white','callback',...
    'SHOWNODES(ShowNodes,ShowElements,revcoordinates,X,Y,r,n,nodes,alpha_field)');
pos = get(ShowNodes,'position') ;
pos = [2*pos(1) pos(2) pos(3) pos(4)] ;
ShowElements = uicontrol('style','toggle','string','Elements','value',0,....
    'position',[pos(1) pos(2) pos(3) pos(4)],'background','white','callback',....
    'SHOWELEMENTS(ShowElements,ShowNodes,revcoordinates,X,Y,iel,nodes,r,alpha_field)');

plot(alpha_field)
hold on


% plot(alpha_inner_field)
% hold on

%{
%To display volume fraction
figure(2)
contourf(nx,ny,T_rev,[1,0]);
xlabel('x[m]','FontSize',16)
ylabel('y[m]','FontSize',16)

title('Fraction Field T','FontSize',16)
grid on
axis equal
colorbar
%}

%bu bölümde bir uyar? var !!!!! cok önemli, iki nokta aras?ndaki dogru
%denklemi her iterasyonda bulunup yazilmali, yoksa yeni gelen de?erler
%eskileri siliyor

    
  % oriented with element number
   inner_nodes_element=zeros(n_active_el,12);
    for nn=1:n_active_el
        
        for n_say=1:6
        dn=nodes(nn,n_say);
   if  ((((revcoordinates(dn,1)-x0)^2)+ ((revcoordinates(dn,2)-y0)^2))^0.5)<radius
     
    inner_nodes_element(nn,(2*n_say)-1)=revcoordinates(dn,1);
    inner_nodes_element(nn,(2*n_say))=revcoordinates(dn,2);
    
 
      end
        end 
        
  
    end
    
    
    
    
    for nn=1:n_active_el
   if inner_nodes_element(nn,1)~=0 & inner_nodes_element(nn,3)~=0 & inner_nodes_element(nn,5)~=0 & inner_nodes_element(nn,7)~=0 & inner_nodes_element(nn,9)~=0 & inner_nodes_element(nn,11)~=0
       
     inner_nodes_element(nn,:)=0;
     
     alpha(nn,1)=1;
       
       
   end
    end
    
    
    
    
    inner_nodes_element_neighbor=inner_nodes_element;
  
    for vi=1:n_active_el
        for vii=1:2:11
            
            
            if vii==1
                vi_in=(vii+1)/2;
                if inner_nodes_element(vi,1)~=0 && (inner_nodes_element(vi,11)==0 | inner_nodes_element(vi,3)==0)
        
                    if inner_nodes_element(vi,11)==0 & inner_nodes_element(vi,3)==0
                        
                                 duu(vi,:) = [revcoordinates(nodes(vi,6),1);revcoordinates(nodes(vi,6),2);inner_nodes_element(vi,vii);inner_nodes_element(vi,vii+1);revcoordinates(nodes(vi,2),1);revcoordinates(nodes(vi,2),2);];    
                                  
    x1=duu(vi,3);
    
    y1=duu(vi,4);
    
    x2=duu(vi,1);
    
    y2=duu(vi,2);
    
    if x1~=x2
    
    m=(y1-y2)/(x1-x2);  % slope of line
        
    
    %this part is comes from equilibrum equations of circle eq. and line
    %eq.  intersection points
    
    
    K=-m*x1+y1-y0;
    
    A=(m^2)+1;
    
    B=-2*x0 +2*m*K;
    
    C=x0^2 +K^2 -radius^2;
    
    delta= B^2 -4*A*C;  %discriminant of equilibrum function
    
    
    rootx1=(-B+(delta^(0.5)))/(2*A);
    
    rootx2=(-B-(delta^(0.5)))/(2*A);
   
   rooty1=m*(rootx1-x1)+y1;
    
    rooty2=m*(rootx2-x1)+y1;
    %arrange with part of circle 
    
    
    r1= (((rootx1-x1)^2) + ((rooty1-y1)^2))^0.5;
    
   r2= (((rootx2-x1)^2) + ((rooty2-y1)^2))^0.5;
    
    if r1<r2
        
        rootx=rootx1;
        rooty=rooty1;
        
    elseif r2<r1
    
        rootx=rootx2;
        rooty=rooty2;
    end
    
                 if inner_nodes_element_neighbor(vi,11)==0

        inner_nodes_element_neighbor(vi,11)=rootx;
        inner_nodes_element_neighbor(vi,12)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,11)~=0
        
        inner_nodes_element_neighbor(vi,11,2)=rootx;
        inner_nodes_element_neighbor(vi,12,2)=rooty;
        
    end
  
    
    elseif x1==x2
        
        
         rootx=x1;
        
        
        rooty1=y0+((radius^2)-(rootx-x0)^2)^0.5;
        rooty2=y0-((radius^2)-(rootx-x0)^2)^0.5;
        
        r_in1=y1-rooty1;
        r_in2=y1-rooty2;
        
        if abs(r_in1)<abs(r_in2)
            
            rooty=rooty1;
            
        elseif abs(r_in2)<abs(r_in1)
        
             rooty=rooty2;
        end
               if inner_nodes_element_neighbor(vi,11)==0

        inner_nodes_element_neighbor(vi,11)=rootx;
        inner_nodes_element_neighbor(vi,12)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,11)~=0
        
        inner_nodes_element_neighbor(vi,11,2)=rootx;
        inner_nodes_element_neighbor(vi,12,2)=rooty;
        
    end
  
      
        
    end
        
        
        
        
        
        
        
    
    x1=duu(vi,3);
    
    y1=duu(vi,4);
    
    x2=duu(vi,5);
    
    y2=duu(vi,6);
    
    if x1~=x2
    
    m=(y1-y2)/(x1-x2);  % slope of line
        
    
    %this part is comes from equilibrum equations of circle eq. and line
    %eq.  intersection points
    
    
    K=-m*x1+y1-y0;
    
    A=(m^2)+1;
    
    B=-2*x0 +2*m*K;
    
    C=x0^2 +K^2 -radius^2;
    
    delta= B^2 -4*A*C;  %discriminant of equlubrium function
    
    
    rootx1=(-B+(delta^(0.5)))/(2*A);
    
    rootx2=(-B-(delta^(0.5)))/(2*A);
   
   rooty1=m*(rootx1-x1)+y1;
    
    rooty2=m*(rootx2-x1)+y1;
    %arrange with part of circle 
    
   r1= (((rootx1-x1)^2) + ((rooty1-y1)^2))^0.5;
    
   r2= (((rootx2-x1)^2) + ((rooty2-y1)^2))^0.5;
    
    if r1<r2
        
        rootx=rootx1;
        rooty=rooty1;
        
    elseif r2<r1
    
        rootx=rootx2;
        rooty=rooty2;
    end
    
          if inner_nodes_element_neighbor(vi,vii+2)==0

        inner_nodes_element_neighbor(vi,vii+2)=rootx;
        inner_nodes_element_neighbor(vi,vii+3)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,vii+2)~=0
        
        inner_nodes_element_neighbor(vi,vii+2,2)=rootx;
        inner_nodes_element_neighbor(vi,vii+3,2)=rooty;
        
    end
  
    
    
    elseif x1==x2
        
        
        rootx=x1;
        
        
        rooty1=y0+((radius^2)-(rootx-x0)^2)^0.5;
        rooty2=y0-((radius^2)-(rootx-x0)^2)^0.5;
        
        r_in1=y1-rooty1;
        r_in2=y1-rooty2;
        
        if abs(r_in1)<abs(r_in2)
            
            rooty=rooty1;
            
        elseif abs(r_in2)<abs(r_in1)
        
             rooty=rooty2;
        end
             if inner_nodes_element_neighbor(vi,vii+2)==0

        inner_nodes_element_neighbor(vi,vii+2)=rootx;
        inner_nodes_element_neighbor(vi,vii+3)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,vii+2)~=0
        
        inner_nodes_element_neighbor(vi,vii+2,2)=rootx;
        inner_nodes_element_neighbor(vi,vii+3,2)=rooty;
        
    end
  
        
    end
    
                                 
                    elseif inner_nodes_element(vi,11)==0
                        
                                duu(vi,:) = [revcoordinates(nodes(vi,6),1);revcoordinates(nodes(vi,6),2);inner_nodes_element(vi,vii);inner_nodes_element(vi,vii+1);0;0];    

    x1=duu(vi,3);
    
    y1=duu(vi,4);
    
    x2=duu(vi,1);
    
    y2=duu(vi,2);
    
    if x1~=x2
        
    
    m=(y1-y2)/(x1-x2);  % slope of line
        
    
    %this part is comes from equilibrum equations of circle eq. and line
    %eq.  intersection points
    
    
    K=-m*x1+y1-y0;
    
    A=(m^2)+1;
    
    B=-2*x0 +2*m*K;
    
    C=x0^2 +K^2 -radius^2;
    
    delta= B^2 -4*A*C;  %discriminant of equlubrium function
    
    
    rootx1=(-B+(delta^(0.5)))/(2*A);
    
    rootx2=(-B-(delta^(0.5)))/(2*A);
   
   rooty1=m*(rootx1-x1)+y1;
    
    rooty2=m*(rootx2-x1)+y1;
    %arrange with part of circle 
    
    
   r1= (((rootx1-x1)^2) + ((rooty1-y1)^2))^0.5;
    
   r2= (((rootx2-x1)^2) + ((rooty2-y1)^2))^0.5;
    
    if r1<r2
        
        rootx=rootx1;
        rooty=rooty1;
        
    elseif r2<r1
    
        rootx=rootx2;
        rooty=rooty2;
    end
    
      

        inner_nodes_element_neighbor(vi,13)=rootx;
        inner_nodes_element_neighbor(vi,14)=rooty;
    
    
  
    
    elseif x1==x2
        
          rootx=x1;
        
        
        rooty1=y0+((radius^2)-(rootx-x0)^2)^0.5;
        rooty2=y0-((radius^2)-(rootx-x0)^2)^0.5;
        
        r_in1=y1-rooty1;
        r_in2=y1-rooty2;
        
        if abs(r_in1)<abs(r_in2)
            
            rooty=rooty1;
            
        elseif abs(r_in2)<abs(r_in1)
        
             rooty=rooty2;
        end
         

        inner_nodes_element_neighbor(vi,13)=rootx;
        inner_nodes_element_neighbor(vi,14)=rooty;
    
    
   
        
        
    
        
     
    end
    
    
                                
                    elseif inner_nodes_element(vi,3)==0
                        
                                duu(vi,:) = [0;0;inner_nodes_element(vi,vii);inner_nodes_element(vi,vii+1);revcoordinates(nodes(vi,2),1);revcoordinates(nodes(vi,2),2);];    

    x1=duu(vi,3);
    
    y1=duu(vi,4);
    
    x2=duu(vi,5);
    
    y2=duu(vi,6);
    
    if x1~=x2
    
    m=(y1-y2)/(x1-x2);  % slope of line
        
    
    %this part is comes from equilibrum equations of circle eq. and line
    %eq.  intersection points
    
    
    K=-m*x1+y1-y0;
    
    A=(m^2)+1;
    
    B=-2*x0 +2*m*K;
    
    C=x0^2 +K^2 -radius^2;
    
    delta= B^2 -4*A*C;  %discriminant of equilibrum function
    
    
    rootx1=(-B+(delta^(0.5)))/(2*A);
    
    rootx2=(-B-(delta^(0.5)))/(2*A);
   
    rooty1=m*(rootx1-x1)+y1;
    
    rooty2=m*(rootx2-x1)+y1;
    %arrange with part of circle 
    
    
   r1= (((rootx1-x1)^2) + ((rooty1-y1)^2))^0.5;
    
   r2= (((rootx2-x1)^2) + ((rooty2-y1)^2))^0.5;
    
    if r1<r2
        
        rootx=rootx1;
        rooty=rooty1;
        
    elseif r2<r1
    
        rootx=rootx2;
        rooty=rooty2;
    end
    
      if inner_nodes_element_neighbor(vi,3)==0

        inner_nodes_element_neighbor(vi,3)=rootx;
        inner_nodes_element_neighbor(vi,4)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,3)~=0
        
        inner_nodes_element_neighbor(vi,3,2)=rootx;
        inner_nodes_element_neighbor(vi,4,2)=rooty;
        
    end
        
    elseif x1==x2
        
          rootx=x1;
        
        
        rooty1=y0+((radius^2)-(rootx-x0)^2)^0.5;
        rooty2=y0-((radius^2)-(rootx-x0)^2)^0.5;
        
        r_in1=y1-rooty1;
        r_in2=y1-rooty2;
        
        if abs(r_in1)<abs(r_in2)
            
            rooty=rooty1;
            
        elseif abs(r_in2)<abs(r_in1)
        
             rooty=rooty2;
        end
         if inner_nodes_element_neighbor(vi,3)==0

        inner_nodes_element_neighbor(vi,3)=rootx;
        inner_nodes_element_neighbor(vi,4)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,3)~=0
        
        inner_nodes_element_neighbor(vi,3,2)=rootx;
        inner_nodes_element_neighbor(vi,4,2)=rooty;
        
    end
        
        
        
    end
                                
                    end
                end
    
                
    elseif vii>=3 & vii<=9
               
    if inner_nodes_element(vi,vii)~=0 & (inner_nodes_element(vi,vii-2)==0 || inner_nodes_element(vi,vii+2)==0)
        vi_in=(vii+1)/2;
        
        %en son buras? yap?ld?  79. eleman?n kom?ular? tamamland?
       if  inner_nodes_element(vi,vii+2)==0 & inner_nodes_element(vi,vii-2)~=0
           
                 duu(vi,:) = [0;0;inner_nodes_element(vi,vii);inner_nodes_element(vi,vii+1);revcoordinates(nodes(vi,vi_in+1),1);revcoordinates(nodes(vi,vi_in+1),2)];    

    x1=duu(vi,3);
    
    y1=duu(vi,4);
    
    x2=duu(vi,5);
    
    y2=duu(vi,6);
    
    if x1~=x2
    
    m=(y1-y2)/(x1-x2);  % slope of line
        
    
    %this part is comes from equilibrum equations of circle eq. and line
    %eq.  intersection points
    
    
    K=-m*x1+y1-y0;
    
    A=(m^2)+1;
    
    B=-2*x0 +2*m*K;
    
    C=x0^2 +K^2 -radius^2;
    
    delta= B^2 -4*A*C;  %discriminant of equlubrium function
    
    
    

    rootx1=(-B+(delta^(0.5)))/(2*A);
    
    rootx2=(-B-(delta^(0.5)))/(2*A);
   
    rooty1=m*(rootx1-x1)+y1;
    
    rooty2=m*(rootx2-x1)+y1;
    %arrange with part of circle 
    
    
   r1= ((abs(rootx1-x1)^2) + (abs(rooty1-y1)^2))^0.5;
    
   r2= ((abs(rootx2-x1)^2) + (abs(rooty2-y1)^2))^0.5;
    
    if r1<r2
        
        rootx=rootx1;
        rooty=rooty1;
        
    elseif r2<r1
    
        rootx=rootx2;
        rooty=rooty2;
        
        
       
    end
    
    if inner_nodes_element_neighbor(vi,vii+2)==0

    inner_nodes_element_neighbor(vi,vii+2)=rootx;
    inner_nodes_element_neighbor(vi,vii+3)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,vii+2)~=0
        
        inner_nodes_element_neighbor(vi,vii+2,2)=rootx;
        inner_nodes_element_neighbor(vi,vii+3,2)=rooty;
        
    end
    
    
    
    elseif x1==x2
        
        rootx=x1;
        
        
        rooty1=y0+((radius^2)-(rootx-x0)^2)^0.5;
        rooty2=y0-((radius^2)-(rootx-x0)^2)^0.5;
        
        r_in1=y1-rooty1;
        r_in2=y1-rooty2;
        
        if abs(r_in1)<abs(r_in2)
            
            rooty=rooty1;
            
        elseif abs(r_in2)<abs(r_in1)
        
             rooty=rooty2;
        end
        
         if inner_nodes_element_neighbor(vi,vii+2)==0

        inner_nodes_element_neighbor(vi,vii+2)=rootx;
        inner_nodes_element_neighbor(vi,vii+3)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,vii+2)~=0
        
        inner_nodes_element_neighbor(vi,vii+2,2)=rootx;
        inner_nodes_element_neighbor(vi,vii+3,2)=rooty;
        
    end
  
        
        
    end
         
           
          
           
       elseif     inner_nodes_element(vi,vii-2)==0 & inner_nodes_element(vi,vii+2)==0 
        
                 duu(vi,:) = [revcoordinates(nodes(vi,vi_in-1),1);revcoordinates(nodes(vi,vi_in-1),2);inner_nodes_element(vi,vii);inner_nodes_element(vi,vii+1);revcoordinates(nodes(vi,vi_in+1),1);revcoordinates(nodes(vi,vi_in+1),2)];
    
    x1=duu(vi,3);
    
    y1=duu(vi,4);
    
    x2=duu(vi,1);
    
    y2=duu(vi,2);
    
    if x1~=x2
        
    
    m=(y1-y2)/(x1-x2);  % slope of line
        
    
    %this part is comes from equilibrum equations of circle eq. and line
    %eq.  intersection points
    
    
    K=-m*x1+y1-y0;
    
    A=(m^2)+1;
    
    B=-2*x0 +2*m*K;
    
    C=x0^2 +K^2 -radius^2;
    
    delta= B^2 -4*A*C;  %discriminant of equilibrum function
    
    
    rootx1=(-B+(delta^(0.5)))/(2*A);
    
    rootx2=(-B-(delta^(0.5)))/(2*A);
   
    rooty1=m*(rootx1-x1)+y1;
    
    rooty2=m*(rootx2-x1)+y1;
    %arrange with part of circle 
    
    
   r1= (((rootx1-x1)^2) + ((rooty1-y1)^2))^0.5;
    
   r2= (((rootx2-x1)^2) + ((rooty2-y1)^2))^0.5;
    
    if r1<r2
        
        rootx=rootx1;
        rooty=rooty1;
        
    elseif r2<r1
    
        rootx=rootx2;
        rooty=rooty2;
    end
    
    if inner_nodes_element_neighbor(vi,vii-2)==0

    inner_nodes_element_neighbor(vi,vii-2)=rootx;
    inner_nodes_element_neighbor(vi,vii-1)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,vii-2)~=0
        
        inner_nodes_element_neighbor(vi,vii-2,2)=rootx;
        inner_nodes_element_neighbor(vi,vii-1,2)=rooty;
        
    end
    
    
    elseif x1==x2
        
          rootx=x1;
        
        
        rooty1=y0+((radius^2)-(rootx-x0)^2)^0.5;
        rooty2=y0-((radius^2)-(rootx-x0)^2)^0.5;
        
        r_in1=y1-rooty1;
        r_in2=y1-rooty2;
        
        if abs(r_in1)<abs(r_in2)
            
            rooty=rooty1;
            
        elseif abs(r_in2)<abs(r_in1)
        
             rooty=rooty2;
        end
        
        if inner_nodes_element_neighbor(vi,vii-2)==0

    inner_nodes_element_neighbor(vi,vii-2)=rootx;
    inner_nodes_element_neighbor(vi,vii-1)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,vii-2)~=0
        
        inner_nodes_element_neighbor(vi,vii-2,2)=rootx;
        inner_nodes_element_neighbor(vi,vii-1,2)=rooty;
        
    end
       
        
    end
        
        

    x1=duu(vi,3);
    
    y1=duu(vi,4);
    
    x2=duu(vi,5);
    
    y2=duu(vi,6);
    
    if x1~=x2
        
    
    m=(y1-y2)/(x1-x2);  % slope of line
        
    
    %this part is comes from equilibrum equations of circle eq. and line
    %eq.  intersection points
    
    
    K=-m*x1+y1-y0;
    
    A=(m^2)+1;
    
    B=-2*x0 +2*m*K;
    
    C=x0^2 +K^2 -radius^2;
    
    delta= B^2 -4*A*C;  %discriminant of equlubrium function
    
    
    rootx1=(-B+(delta^(0.5)))/(2*A);
    
    rootx2=(-B-(delta^(0.5)))/(2*A);
   
   rooty1=m*(rootx1-x1)+y1;
    
    rooty2=m*(rootx2-x1)+y1;
    %arrange with part of circle 
    
   r1= (((rootx1-x1)^2) + ((rooty1-y1)^2))^0.5;
    
   r2= (((rootx2-x1)^2) + ((rooty2-y1)^2))^0.5;
    
    if r1<r2
        
        rootx=rootx1;
        rooty=rooty1;
        
    elseif r2<r1
    
        rootx=rootx2;
        rooty=rooty2;
    end
    if inner_nodes_element_neighbor(vi,vii+2)==0

    inner_nodes_element_neighbor(vi,vii+2)=rootx;
    inner_nodes_element_neighbor(vi,vii+3)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,vii+2)~=0
        
        inner_nodes_element_neighbor(vi,vii+2,2)=rootx;
        inner_nodes_element_neighbor(vi,vii+3,2)=rooty;
        
    end
    
      
    
    elseif x1==x2
        
          rootx=x1;
        
        
        rooty1=y0+((radius^2)-(rootx-x0)^2)^0.5;
        rooty2=y0-((radius^2)-(rootx-x0)^2)^0.5;
        
        r_in1=y1-rooty1;
        r_in2=y1-rooty2;
        
        if abs(r_in1)<abs(r_in2)
            
            rooty=rooty1;
            
        elseif abs(r_in2)<abs(r_in1)
        
             rooty=rooty2;
        end
        
        if inner_nodes_element_neighbor(vi,vii+2)==0

    inner_nodes_element_neighbor(vi,vii+2)=rootx;
    inner_nodes_element_neighbor(vi,vii+3)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,vii+2)~=0
        
        inner_nodes_element_neighbor(vi,vii+2,2)=rootx;
        inner_nodes_element_neighbor(vi,vii+3,2)=rooty;
        
        end
            
    end
    
   

       elseif    inner_nodes_element(vi,vii-2)==0 
                    
                 duu(vi,:) = [revcoordinates(nodes(vi,vi_in-1),1);revcoordinates(nodes(vi,vi_in-1),2);inner_nodes_element(vi,vii);inner_nodes_element(vi,vii+1);0;0];    
    
    x1=duu(vi,3);
    
    y1=duu(vi,4);
    
    x2=duu(vi,1);
    
    y2=duu(vi,2);
    
    if x1~=x2
       
    
    m=(y1-y2)/(x1-x2);  % slope of line
        
    
    %this part is comes from equilibrum equations of circle eq. and line
    %eq.  intersection points
    
    
    K=-m*x1+y1-y0;
    
    A=(m^2)+1;
    
    B=-2*x0 +2*m*K;
    
    C=x0^2 +K^2 -radius^2;
    
    delta= B^2 -4*A*C;  %discriminant of equlubrium function
    
    
    rootx1=(-B+(delta^(0.5)))/(2*A);
    
    rootx2=(-B-(delta^(0.5)))/(2*A);
   
   rooty1=m*(rootx1-x1)+y1;
    
    rooty2=m*(rootx2-x1)+y1;
    %arrange with part of circle 
    
    
  r1= (((rootx1-x1)^2) + ((rooty1-y1)^2))^0.5;
    
   r2= (((rootx2-x1)^2) + ((rooty2-y1)^2))^0.5;
    
    if r1<r2
        
        rootx=rootx1;
        rooty=rooty1;
        
    elseif r2<r1
    
        rootx=rootx2;
        rooty=rooty2;
    end
    
    if inner_nodes_element_neighbor(vi,vii-2)==0

        inner_nodes_element_neighbor(vi,vii-2)=rootx;
        inner_nodes_element_neighbor(vi,vii-1)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,vii-2)~=0
        
        inner_nodes_element_neighbor(vi,vii-2,2)=rootx;
        inner_nodes_element_neighbor(vi,vii-1,2)=rooty;
        
    end
    
    
    elseif x1==x2
        
          rootx=x1;
        
        
        rooty1=y0+((radius^2)-(rootx-x0)^2)^0.5;
        rooty2=y0-((radius^2)-(rootx-x0)^2)^0.5;
        
        r_in1=y1-rooty1;
        r_in2=y1-rooty2;
        
        if abs(r_in1)<abs(r_in2)
            
            rooty=rooty1;
            
        elseif abs(r_in2)<abs(r_in1)
        
             rooty=rooty2;
        end
        
        if inner_nodes_element_neighbor(vi,vii-2)==0

        inner_nodes_element_neighbor(vi,vii-2)=rootx;
        inner_nodes_element_neighbor(vi,vii-1)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,vii-2)~=0
        
        inner_nodes_element_neighbor(vi,vii-2,2)=rootx;
        inner_nodes_element_neighbor(vi,vii-1,2)=rooty;
        
        end
       
  
           
    end
    
    
%     if vi==95 && vii==7
%         rootx_sec=[rootx;rooty;x1;y1;x2;y2;r1;r2;rootx1;rooty1;rootx2;rooty2;A;B;C;delta;K;m]
%        
%     end
                 
                 
          
       end
       
    end
           
                
                
            elseif vii==11
                vi_in=(vii+1)/2;
                if inner_nodes_element(vi,11)~=0 && (inner_nodes_element(vi,9)==0 | inner_nodes_element(vi,1)==0)
                    
                    if inner_nodes_element(vi,9)==0 & inner_nodes_element(vi,1)==0
        
                    duu(vi,:) = [revcoordinates(nodes(vi,5),1);revcoordinates(nodes(vi,5),2);inner_nodes_element(vi,vii);inner_nodes_element(vi,vii+1);revcoordinates(nodes(vi,1),1);revcoordinates(nodes(vi,1),2);];    
   
                    
    x1=duu(vi,3);
    
    y1=duu(vi,4);
    
    x2=duu(vi,5);
    
    y2=duu(vi,6);
    
    if x1~=x2
    
    m=(y1-y2)/(x1-x2);  % slope of line
        
    
    %this part is comes from equilibrum equations of circle eq. and line
    %eq.  intersection points
    
    
    K=-m*x1+y1-y0;
    
    A=(m^2)+1;
    
    B=-2*x0 +2*m*K;
    
    C=x0^2 +K^2 -radius^2;
    
    delta= B^2 -4*A*C;  %discriminant of equlubrium function
    
    
    rootx1=(-B+(delta^(0.5)))/(2*A);
    
    rootx2=(-B-(delta^(0.5)))/(2*A);
   
   rooty1=m*(rootx1-x1)+y1;
    
    rooty2=m*(rootx2-x1)+y1;
    %arrange with part of circle 
    
    
   r1= (((rootx1-x1)^2) + ((rooty1-y1)^2))^0.5;
    
   r2= (((rootx2-x1)^2) + ((rooty2-y1)^2))^0.5;
    
    if r1<r2
        
        rootx=rootx1;
        rooty=rooty1;
        
    elseif r2<r1
    
        rootx=rootx2;
        rooty=rooty2;
    end
    
     if inner_nodes_element_neighbor(vi,1)==0

        inner_nodes_element_neighbor(vi,1)=rootx;
        inner_nodes_element_neighbor(vi,2)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,1)~=0
        
        inner_nodes_element_neighbor(vi,1,2)=rootx;
        inner_nodes_element_neighbor(vi,2,2)=rooty;
        
    end 
    elseif x1==x2
        
          rootx=x1;
        
        
        rooty1=y0+((radius^2)-(rootx-x0)^2)^0.5;
        rooty2=y0-((radius^2)-(rootx-x0)^2)^0.5;
        
        r_in1=y1-rooty1;
        r_in2=y1-rooty2;
        
        if abs(r_in1)<abs(r_in2)
            
            rooty=rooty1;
            
        elseif abs(r_in2)<abs(r_in1)
        
             rooty=rooty2;
        end
        if inner_nodes_element_neighbor(vi,1)==0

        inner_nodes_element_neighbor(vi,1)=rootx;
        inner_nodes_element_neighbor(vi,2)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,1)~=0
        
        inner_nodes_element_neighbor(vi,1,2)=rootx;
        inner_nodes_element_neighbor(vi,2,2)=rooty;
        
    end
      
        
    end
                    
    x1=duu(vi,3);
    
    y1=duu(vi,4);
    
    x2=duu(vi,1);
    
    y2=duu(vi,2);
    
    if x1~=x2
        
        
    m=(y1-y2)/(x1-x2);  % slope of line
        
    
    %this part is comes from equilibrum equations of circle eq. and line
    %eq.  intersection points
    
    
    K=-m*x1+y1-y0;
    
    A=(m^2)+1;
    
    B=-2*x0 +2*m*K;
    
    C=x0^2 +K^2 -radius^2;
    
    delta= B^2 -4*A*C;  %discriminant of equlubrium function
    
    
    rootx1=(-B+(delta^(0.5)))/(2*A);
    
    rootx2=(-B-(delta^(0.5)))/(2*A);
   
   rooty1=m*(rootx1-x1)+y1;
    
    rooty2=m*(rootx2-x1)+y1;
    %arrange with part of circle 
    
    
  r1= (((rootx1-x1)^2) + ((rooty1-y1)^2))^0.5;
    
   r2= (((rootx2-x1)^2) + ((rooty2-y1)^2))^0.5;
    
    if r1<r2
        
        rootx=rootx1;
        rooty=rooty1;
        
    elseif r2<r1
    
        rootx=rootx2;
        rooty=rooty2;
    end
   if inner_nodes_element_neighbor(vi,vii-2)==0

        inner_nodes_element_neighbor(vi,vii-2)=rootx;
        inner_nodes_element_neighbor(vi,vii-1)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,vii-2)~=0
        
        inner_nodes_element_neighbor(vi,vii-2,2)=rootx;
        inner_nodes_element_neighbor(vi,vii-1,2)=rooty;
        
    end
                    
    elseif x1==x2
          rootx=x1;
        
        
        rooty1=y0+((radius^2)-(rootx-x0)^2)^0.5;
        rooty2=y0-((radius^2)-(rootx-x0)^2)^0.5;
        
        r_in1=y1-rooty1;
        r_in2=y1-rooty2;
        
        if abs(r_in1)<abs(r_in2)
            
            rooty=rooty1;
            
        elseif abs(r_in2)<abs(r_in1)
        
             rooty=rooty2;
        end
        if inner_nodes_element_neighbor(vi,vii-2)==0

        inner_nodes_element_neighbor(vi,vii-2)=rootx;
        inner_nodes_element_neighbor(vi,vii-1)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,vii-2)~=0
        
        inner_nodes_element_neighbor(vi,vii-2,2)=rootx;
        inner_nodes_element_neighbor(vi,vii-1,2)=rooty;
        
    end
    
      
    end
                    
                    elseif inner_nodes_element(vi,9)==0
                        
                            duu(vi,:) = [revcoordinates(nodes(vi,5),1);revcoordinates(nodes(vi,5),2);inner_nodes_element(vi,vii);inner_nodes_element(vi,vii+1);0;0];    

                    
    x1=duu(vi,3);
    
    y1=duu(vi,4);
    
    x2=duu(vi,1);
    
    y2=duu(vi,2);
    
    if x1~=x2
        
    
    m=(y1-y2)/(x1-x2);  % slope of line
        
    
    %this part is comes from equilibrum equations of circle eq. and line
    %eq.  intersection points
    
    
    K=-m*x1+y1-y0;
    
    A=(m^2)+1;
    
    B=-2*x0 +2*m*K;
    
    C=x0^2 +K^2 -radius^2;
    
    delta= B^2 -4*A*C;  %discriminant of equlubrium function
    
    
    rootx1=(-B+(delta^(0.5)))/(2*A);
    
    rootx2=(-B-(delta^(0.5)))/(2*A);
   
   rooty1=m*(rootx1-x1)+y1;
    
    rooty2=m*(rootx2-x1)+y1;
    %arrange with part of circle 
    
    
   r1= (((rootx1-x1)^2) + ((rooty1-y1)^2))^0.5;
    
   r2= (((rootx2-x1)^2) + ((rooty2-y1)^2))^0.5;
    
    if r1<r2
        
        rootx=rootx1;
        rooty=rooty1;
        
    elseif r2<r1
    
        rootx=rootx2;
        rooty=rooty2;
    end
    
      if inner_nodes_element_neighbor(vi,vii-2)==0

        inner_nodes_element_neighbor(vi,vii-2)=rootx;
        inner_nodes_element_neighbor(vi,vii-1)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,vii-2)~=0
        
        inner_nodes_element_neighbor(vi,vii-2,2)=rootx;
        inner_nodes_element_neighbor(vi,vii-1,2)=rooty;
        
    end
    

    elseif x1==x2
        
          rootx=x1;
        
        
        rooty1=y0+((radius^2)-(rootx-x0)^2)^0.5;
        rooty2=y0-((radius^2)-(rootx-x0)^2)^0.5;
        
        r_in1=y1-rooty1;
        r_in2=y1-rooty2;
        
        if abs(r_in1)<abs(r_in2)
            
            rooty=rooty1;
            
        elseif abs(r_in2)<abs(r_in1)
        
             rooty=rooty2;
        end
        
          if inner_nodes_element_neighbor(vi,vii-2)==0

        inner_nodes_element_neighbor(vi,vii-2)=rootx;
        inner_nodes_element_neighbor(vi,vii-1)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,vii-2)~=0
        
        inner_nodes_element_neighbor(vi,vii-2,2)=rootx;
        inner_nodes_element_neighbor(vi,vii-1,2)=rooty;
        
    end
  
        
    end
                            
                            
                            
                    elseif inner_nodes_element(vi,1)==0
                    
                                    duu(vi,:) = [0;0;inner_nodes_element(vi,vii);inner_nodes_element(vi,vii+1);revcoordinates(nodes(vi,1),1);revcoordinates(nodes(vi,1),2);];

                    
    x1=duu(vi,3);
    
    y1=duu(vi,4);
    
    x2=duu(vi,5);
    
    y2=duu(vi,6);
    
    if x1~=x2
       
    
    m=(y1-y2)/(x1-x2);  % slope of line
        
    
    %this part is comes from equilibrum equations of circle eq. and line
    %eq.  intersection points
    
    
    K=-m*x1+y1-y0;
    
    A=(m^2)+1;
    
    B=-2*x0 +2*m*K;
    
    C=x0^2 +K^2 -radius^2;
    
    delta= B^2 -4*A*C;  %discriminant of equlubrium function
    
    
    rootx1=(-B+(delta^(0.5)))/(2*A);
    
    rootx2=(-B-(delta^(0.5)))/(2*A);
   
   rooty1=m*(rootx1-x1)+y1;
    
    rooty2=m*(rootx2-x1)+y1;
    %arrange with part of circle 
    
    
  r1= (((rootx1-x1)^2) + ((rooty1-y1)^2))^0.5;
    
   r2= (((rootx2-x1)^2) + ((rooty2-y1)^2))^0.5;
    
    if r1<r2
        
        rootx=rootx1;
        rooty=rooty1;
        
    elseif r2<r1
    
        rootx=rootx2;
        rooty=rooty2;
    end
    
      if inner_nodes_element_neighbor(vi,1)==0

        inner_nodes_element_neighbor(vi,1)=rootx;
        inner_nodes_element_neighbor(vi,2)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,1)~=0
        
        inner_nodes_element_neighbor(vi,1,2)=rootx;
        inner_nodes_element_neighbor(vi,2,2)=rooty;
        
    end
    
   
                                    
    elseif x1==x2
        
          rootx=x1;
        
        
        rooty1=y0+((radius^2)-(rootx-x0)^2)^0.5;
        rooty2=y0-((radius^2)-(rootx-x0)^2)^0.5;
        
        r_in1=y1-rooty1;
        r_in2=y1-rooty2;
        
        if abs(r_in1)<abs(r_in2)
            
            rooty=rooty1;
            
        elseif abs(r_in2)<abs(r_in1)
        
             rooty=rooty2;
        end
        if inner_nodes_element_neighbor(vi,1)==0

        inner_nodes_element_neighbor(vi,1)=rootx;
        inner_nodes_element_neighbor(vi,2)=rooty;
    
    
    elseif inner_nodes_element_neighbor(vi,1)~=0
        
        inner_nodes_element_neighbor(vi,1,2)=rootx;
        inner_nodes_element_neighbor(vi,2,2)=rooty;
        
    end
    
        
        
    end
    
                    end
    
    end
            end
                
                
        end
        
    end
    
 
    
     
   
    
    matrix1=inner_nodes_element_neighbor(:,:,1);
    matrix2=inner_nodes_element_neighbor(:,:,2);
   
    
    
    
       
   
    
    
    
    for is=1:iel
    
    side = nnz(inner_nodes_element_neighbor(is,:,1));
    side2= nnz(inner_nodes_element_neighbor(is,:,2));
    
    if side2~=0
       side_n=(side+side2)/2; 
       
    else
            side_n=side/2;
        
    end
    

    
    if side_n==7
        
        scd=1;
        scdn=1;
        for sc=1:6
        
            
       if inner_nodes_element_neighbor(is,sc*2-1,1)~=0
           
           if sc==1 & inner_nodes_element_neighbor(is,13,1)~=0
                
                P(scd,1)=inner_nodes_element_neighbor(is,13,1);
                P(scd,2)=inner_nodes_element_neighbor(is,14,1);
                P(scd+1,1)=inner_nodes_element_neighbor(is,2*scdn-1,1);
                P(scd+1,2)=inner_nodes_element_neighbor(is,2*scdn,1);
                
                
                scd=scd+2;
                scdn=scdn+1;
                
           elseif sc==1 & inner_nodes_element_neighbor(is,1,2)~=0
               
                P(scd,1)=inner_nodes_element_neighbor(is,2*scdn-1,2);
                P(scd,2)=inner_nodes_element_neighbor(is,2*scdn,2);
                P(scd+1,1)=inner_nodes_element_neighbor(is,2*scdn-1,1);
                P(scd+1,2)=inner_nodes_element_neighbor(is,2*scdn,1);
                
                
                scd=scd+2;
                scdn=scdn+1;
                
           else
               
               if inner_nodes_element_neighbor(is,sc*2-1,2)~=0
                   
               
                    P(scd,1)=inner_nodes_element_neighbor(is,2*scdn-1,1);
                    P(scd,2)=inner_nodes_element_neighbor(is,2*scdn,1);
                    P(scd+1,1)=inner_nodes_element_neighbor(is,2*scdn-1,2);
                    P(scd+1,2)=inner_nodes_element_neighbor(is,2*scdn,2);
               
                    scd=scd+2;
                    scdn=scdn+1;
                    
               elseif  inner_nodes_element_neighbor(is,sc*2-1,2)==0
                    P(scd,1)=inner_nodes_element_neighbor(is,2*scdn-1,1);
                    P(scd,2)=inner_nodes_element_neighbor(is,2*scdn,1);
               
                    scd=scd+1;
                    scdn=scdn+1;
                    
               end
                   
                   
           end
       end
       
            
        end
        
            
   
            polyin(is)= polyshape(P);
            alpha(is,1)=polyin(is).area/general_area;
            
            
            
        
            elseif side_n==6    
                
                P=zeros(6,2);
                
        scd=1;
        scdn=1;
        for sc=1:6
                                  
                    if inner_nodes_element_neighbor(is,sc*2-1,1)~=0 && (scd <=side_n)
                        
                        
                     if sc==1 & inner_nodes_element_neighbor(is,13,1)~=0
                
                    P(scd,1)=inner_nodes_element_neighbor(is,13,1);
                    P(scd,2)=inner_nodes_element_neighbor(is,14,1);
                    P(scd+1,1)=inner_nodes_element_neighbor(is,2*scdn-1,1);
                    P(scd+1,2)=inner_nodes_element_neighbor(is,2*scdn,1);
                
                
                    scd=scd+2;
                    scdn=scdn+1;
                
                    else
                            
                    P(scd,1)=inner_nodes_element_neighbor(is,2*scdn-1,1);
                    P(scd,2)=inner_nodes_element_neighbor(is,2*scdn,1);
               
                    scd=scd+1;
                    scdn=scdn+1;
                       
                     end  
                     
                    else
                         scdn=scdn+1;
                    end
        end

         polyin(is)= polyshape(P);
         alpha(is,1)=polyin(is).area/general_area;

    elseif side_n==5
        
         P=zeros(5,2);
                        
           scd=1;
        scdn=1;
        for sc=1:6
                                  
                    if inner_nodes_element_neighbor(is,sc*2-1,1)~=0 && (scd <=side_n)
                        
                      if sc==1 & inner_nodes_element_neighbor(is,13,1)~=0
                
                    P(scd,1)=inner_nodes_element_neighbor(is,13,1);
                    P(scd,2)=inner_nodes_element_neighbor(is,14,1);
                    P(scd+1,1)=inner_nodes_element_neighbor(is,2*scdn-1,1);
                    P(scd+1,2)=inner_nodes_element_neighbor(is,2*scdn,1);
                
                
                    scd=scd+2;
                    scdn=scdn+1;
                
                    else
                            
                    P(scd,1)=inner_nodes_element_neighbor(is,2*scdn-1,1);
                    P(scd,2)=inner_nodes_element_neighbor(is,2*scdn,1);
               
                    scd=scd+1;
                    scdn=scdn+1;
                    
                                
                    end          
                    else
                        scdn=scdn+1;
                        
                    end
        end

         polyin(is)= polyshape(P);
         alpha(is,1)=polyin(is).area/general_area;
        
        
        
    elseif side_n==4
        
         P=zeros(4,2);
         
           scd=1;
        scdn=1;
        for sc=1:6
                                  
                    if inner_nodes_element_neighbor(is,sc*2-1,1)~=0 && (scd <=side_n)
                        
                      if sc==1 & inner_nodes_element_neighbor(is,13,1)~=0
                
                    P(scd,1)=inner_nodes_element_neighbor(is,13,1);
                    P(scd,2)=inner_nodes_element_neighbor(is,14,1);
                    P(scd+1,1)=inner_nodes_element_neighbor(is,2*scdn-1,1);
                    P(scd+1,2)=inner_nodes_element_neighbor(is,2*scdn,1);
                
                
                    scd=scd+2;
                    scdn=scdn+1;
                
                    else
                            
                    P(scd,1)=inner_nodes_element_neighbor(is,2*scdn-1,1);
                    P(scd,2)=inner_nodes_element_neighbor(is,2*scdn,1);
               
                    scd=scd+1;
                    scdn=scdn+1;
                       
                      end 
                     
                       else
                        scdn=scdn+1;
                                      
                    end
        end

         polyin(is)= polyshape(P);
         alpha(is,1)=polyin(is).area/general_area;
        
    elseif side_n==3
        
         P=zeros(3,2);
           scd=1;
        scdn=1;
        for sc=1:6
                                  
                    if inner_nodes_element_neighbor(is,sc*2-1,1)~=0 && (scd <=side_n)
                        
                      if sc==1 & inner_nodes_element_neighbor(is,13,1)~=0
                
                    P(scd,1)=inner_nodes_element_neighbor(is,13,1);
                    P(scd,2)=inner_nodes_element_neighbor(is,14,1);
                    P(scd+1,1)=inner_nodes_element_neighbor(is,2*scdn-1,1);
                    P(scd+1,2)=inner_nodes_element_neighbor(is,2*scdn,1);
                
                
                    scd=scd+2;
                    scdn=scdn+1;
                
                    else
                            
                    P(scd,1)=inner_nodes_element_neighbor(is,2*scdn-1,1);
                    P(scd,2)=inner_nodes_element_neighbor(is,2*scdn,1);
               
                    scd=scd+1;
                    scdn=scdn+1;
                       
                      end  
                     
                      else
                        scdn=scdn+1;
                                      
                    end
        end

         polyin(is)= polyshape(P);
   
         alpha(is,1)=polyin(is).area/general_area;
    
        
    end
        
    
    end
    
    
    
    
   
    
 %Distributing all variable //node number // surface normal ..

 dist_snormal_alpha=[dist(:,1) snormal(:,1:2) s_value(:,1) dist(:,2) snormal(:,3:4) s_value(:,2) dist(:,3) snormal(:,5:6) s_value(:,3) dist(:,4) snormal(:,7:8) s_value(:,4) dist(:,5) snormal(:,9:10) s_value(:,5) dist(:,6) snormal(:,11:12) s_value(:,6) dist(:,1) alpha  ];
 
 
 
     %surface value calculation
 for el = lim_l:lim_u
  
     
    check1=ismember(el,nbos_el);
  
    if (check1==1)
        
    else
        
   %oriented by counter-clockwise direction      
   mid=el;
   mid_sol=el-1;
   mid_sol_ust=el+Nx-1;
   mid_sag_ust=el+Nx;
   mid_sag=el+1;
   mid_sag_alt=el-Nx+1  ;
   mid_sol_alt=el-Nx   ;
   
s_value(mid,1)=(dist_snormal_alpha(mid,26)+dist_snormal_alpha(mid_sol,26))/2;
s_value(mid,2)=(dist_snormal_alpha(mid,26)+dist_snormal_alpha(mid_sol_ust,26))/2;
s_value(mid,3)=(dist_snormal_alpha(mid,26)+dist_snormal_alpha(mid_sag_ust,26))/2;
s_value(mid,4)=(dist_snormal_alpha(mid,26)+dist_snormal_alpha(mid_sag,26))/2;   
s_value(mid,5)=(dist_snormal_alpha(mid,26)+dist_snormal_alpha(mid_sag_alt,26))/2;
s_value(mid,6)=(dist_snormal_alpha(mid,26)+dist_snormal_alpha(mid_sol_alt,26))/2;        
    end
       
     

 end

   dist_snormal_alpha=[dist(:,1) snormal(:,1:2) s_value(:,1) dist(:,2) snormal(:,3:4) s_value(:,2) dist(:,3) snormal(:,5:6) s_value(:,3) dist(:,4) snormal(:,7:8) s_value(:,4) dist(:,5) snormal(:,9:10) s_value(:,5) dist(:,6) snormal(:,11:12) s_value(:,6) dist(:,1) alpha  ];

 
 
 %Gradient calculation 

 %attention! there are no area and volume value, should add these variable
 
 for n = 1:(Nx^2)-(Nx/2)
      
  gradient_of_element(n,1)= ((dist_snormal_alpha(n,4)*dist_snormal_alpha(n,2)) + (dist_snormal_alpha(n,8)*dist_snormal_alpha(n,6))+ (dist_snormal_alpha(n,12)*dist_snormal_alpha(n,10))+ (dist_snormal_alpha(n,16)*dist_snormal_alpha(n,14))+ (dist_snormal_alpha(n,20)*dist_snormal_alpha(n,18))+ (dist_snormal_alpha(n,24)*dist_snormal_alpha(n,22)));
 
  
  gradient_of_element(n,2)= ((dist_snormal_alpha(n,4)*dist_snormal_alpha(n,3)) + (dist_snormal_alpha(n,8)*dist_snormal_alpha(n,7))+ (dist_snormal_alpha(n,12)*dist_snormal_alpha(n,11))+ (dist_snormal_alpha(n,16)*dist_snormal_alpha(n,15))+ (dist_snormal_alpha(n,20)*dist_snormal_alpha(n,19))+ (dist_snormal_alpha(n,24)*dist_snormal_alpha(n,23)));

  gradient_of_element(n,3)= (gradient_of_element(n,1)^2 +gradient_of_element(n,2)^2)^0.5 ;

 end
 
 
  gradient_of_element_nondim(:,1)= gradient_of_element(:,1)*(1/l_unit);
 
  
  gradient_of_element_nondim(:,2)= gradient_of_element(:,2)*(1/l_unit);
  
  for n = 1:(Nx^2)-(Nx/2)

  gradient_of_element_nondim(n,3)= (gradient_of_element_nondim(n,1)^2 +gradient_of_element_nondim(n,2)^2)^0.5 ;
 
  end
  %find the local centeral coordinates and put grad matrix 4 and 5 columns
  for n_el=1:iel
       
     sec_g(1,1)=dist(n_el,1);
     sec_g(2,1)=dist(n_el,2);
     sec_g(3,1)=dist(n_el,3);
     sec_g(4,1)=dist(n_el,4);
     sec_g(5,1)=dist(n_el,5);
     sec_g(6,1)=dist(n_el,6);
     
     gradient_of_element(n_el,4)=(revcoordinates(sec_g(1,1),1)+revcoordinates(sec_g(2,1),1)+revcoordinates(sec_g(3,1),1)+revcoordinates(sec_g(4,1),1)+revcoordinates(sec_g(5,1),1)+revcoordinates(sec_g(6,1),1))/6;
     gradient_of_element(n_el,5)=(revcoordinates(sec_g(1,1),2)+revcoordinates(sec_g(2,1),2)+revcoordinates(sec_g(3,1),2)+revcoordinates(sec_g(4,1),2)+revcoordinates(sec_g(5,1),2)+revcoordinates(sec_g(6,1),2))/6;
      
       
     gradient_of_element_nondim(n_el,4)=(revcoordinates(sec_g(1,1),1)+revcoordinates(sec_g(2,1),1)+revcoordinates(sec_g(3,1),1)+revcoordinates(sec_g(4,1),1)+revcoordinates(sec_g(5,1),1)+revcoordinates(sec_g(6,1),1))/6;
     gradient_of_element_nondim(n_el,5)=(revcoordinates(sec_g(1,1),2)+revcoordinates(sec_g(2,1),2)+revcoordinates(sec_g(3,1),2)+revcoordinates(sec_g(4,1),2)+revcoordinates(sec_g(5,1),2)+revcoordinates(sec_g(6,1),2))/6;
      
      
           
       
  end
 
  
%find max value of gradient vector
[row, col] = find(ismember(round(gradient_of_element(:,3),3), max(round(gradient_of_element(:,3),3))));


plot(polyin);
axis equal


figure
plot(alpha_field)
hold on

quiverwcolorbar(gradient_of_element(:,4),gradient_of_element(:,5),gradient_of_element(:,1),gradient_of_element(:,2),0.15,jet) % scale factor must be 0-0.2 range

% colormap jet
% colorbar
title (['Gradient Distribution - Unit Lenght Of Elements:',num2str(l_unit)]); 
x0=700;
y0=300;
width=625;
height=500
set(gcf,'position',[x0,y0,width,height],'FontSize',16)
plot(polyin);
axis equal




% figure
% plot(alpha_field)
% hold on
% quiverwcolorbar(gradient_of_element_nondim(:,4),gradient_of_element_nondim(:,5),gradient_of_element_nondim(:,1),gradient_of_element_nondim(:,2),0.005,jet) % scale factor must be 0-0.2 range
% colormap jet
% colorbar
% title ('Gradient Distribution')
%     
% plot(polyin);
% axis equal



figure
plot(alpha_field)
hold on
quiver(gradient_of_element(row(:),4),gradient_of_element(row(:),5),gradient_of_element(row(:),1),gradient_of_element(row(:),2),0.5,'r')
title ('Maximum Gradient Vector Position')



for n=1:iel
    
   
    if gradient_of_element(n,4)>=x0 && gradient_of_element(n,5)>=y0
    theta(n,1)=atand(abs(gradient_of_element(n,5)-y0)/(abs(gradient_of_element(n,4)-x0)));
    
    elseif gradient_of_element(n,4)<x0 && gradient_of_element(n,5)>y0
        
        theta(n,1)=180-atand(abs(gradient_of_element(n,5)-y0)/(abs(gradient_of_element(n,4)-x0)));
        
    elseif gradient_of_element(n,4)<x0 && gradient_of_element(n,5)<y0
        
                theta(n,1)=180+atand(abs(gradient_of_element(n,5)-y0)/(abs(gradient_of_element(n,4)-x0)));
                
                
    elseif gradient_of_element(n,4)>x0 && gradient_of_element(n,5)<y0
        
            theta(n,1)=360-atand(abs(gradient_of_element(n,5)-y0)/(abs(gradient_of_element(n,4)-x0)));

    end

%     theta(n,1)=atand((abs(gradient_of_element(n,5))-abs(y0))/(abs(gradient_of_element(n,4))-abs(x0)));
        
    r_dist(n,1)=((gradient_of_element(n,5)-y0)^2 + (gradient_of_element(n,4)-x0)^2)^0.5;
    
end


gradient_of_element(:,6)=theta(:,1);  % theta degree of each element

gradient_of_element(:,7)=r_dist(:,1); % centeral distance of each element

figure
scatter(gradient_of_element(:,6),gradient_of_element(:,3),'filled');
xlabel('Theta Degree')
ylabel('Magnitude Of Gradient Vector')

%normalized to gradient value according to max value


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













