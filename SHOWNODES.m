function SHOWNODES(ShowNodes,ShowElements,coordinates,X,Y,nnode,nel,nodes,alpha_field) 
%--------------------------------------------------------------------------
%   Purpose:
%           To display and undisplay the node numbers
%--------------------------------------------------------------------------
% To display Node Numbers
ShowNodes = get(ShowNodes,'Value') ;
ShowElements = get(ShowElements,'Value') ;
% Display only Node Numbers
if ShowNodes==1 && ShowElements==0 
    cla(gcf)
    patch(X,Y,'w') ;
    for i = 1:nnode
        text(coordinates(i,1),coordinates(i,2),int2str(i),'fontsize',8,....
            'fontweight','bold','Color','r'); 
    end  
    plot(alpha_field)
  hold on
end
% Display only Element Numbers
if ShowNodes==0 && ShowElements==1
    cla(gcf) ;
    patch(X,Y,'w') ;
    for i = 1:nel
        EX = coordinates(nodes(i,:),1) ;EY = coordinates(nodes(i,:),2) ;
        pos = [sum(EX)/6,sum(EY)/6] ;
        text(pos(1),pos(2),int2str(i),'fontsize',8, 'fontweight','bold','color','b');
    end
    plot(alpha_field)
  hold on
end
% Display both Node and Element Numbers
if ShowNodes==1 && ShowElements==1
    cla(gcf) ;
    patch(X,Y,'w') ;
    for i = 1:nnode
        text(coordinates(i,1),coordinates(i,2),int2str(i),'fontsize',8,....
            'fontweight','bold','Color','r'); 
    end  
    for i = 1:nel
        EX = coordinates(nodes(i,:),1) ;EY = coordinates(nodes(i,:),2) ;
        pos = [sum(EX)/6,sum(EY)/6] ;
        text(pos(1),pos(2),int2str(i),'fontsize',8, 'fontweight','bold','color','b');
    end
    plot(alpha_field)
  hold on
end
% Display None
if ShowNodes==0 && ShowElements==0
    cla(gcf) ;
    patch(X,Y,'w') ;
    plot(alpha_field)
  hold on
end
    
return