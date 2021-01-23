function SHOWELEMENTS(ShowElements,ShowNodes,revcoordinates,X,Y,iel,nodes,r,alpha_field)
%--------------------------------------------------------------------------
%   Purpose:
%           To display and undisplay the Element numbers
%--------------------------------------------------------------------------
ShowElements = get(ShowElements,'Value') ;
ShowNodes = get(ShowNodes,'Value') ;  
% Diaply only Element Numbers
if ShowElements==1 && ShowNodes==0  
    cla(gcf) ;
    patch(X,Y,'w') ;
    for i = 1:iel
        EX = revcoordinates(nodes(i,:),1) ;EY = revcoordinates(nodes(i,:),2) ;
        pos = [sum(EX)/6,sum(EY)/6] ;
        text(pos(1),pos(2),int2str(i),'fontsize',8, 'fontweight','bold','color','b');
    end
  plot(alpha_field)
  hold on
end
% Display only Node Numbers
if ShowElements==0 && ShowNodes==1
    cla(gcf) ;
    patch(X,Y,'w')
    for i = 1:r
        text(revcoordinates(i,1),revcoordinates(i,2),int2str(i),'fontsize',8,....
            'fontweight','bold','Color','r'); 
    end  
    plot(alpha_field)
  hold on
end
% Display both Element Numbers and Node Numbers
if ShowElements==1 && ShowNodes==1
    cla(gcf) ;
    patch(X,Y,'w') ;
    for i = 1:iel
        EX = revcoordinates(nodes(i,:),1) ;EY = revcoordinates(nodes(i,:),2) ;
        pos = [sum(EX)/6,sum(EY)/6] ;
        text(pos(1),pos(2),int2str(i),'fontsize',8, 'fontweight','bold','color','b');
    end
    for i = 1:r
        text(revcoordinates(i,1),revcoordinates(i,2),int2str(i),'fontsize',8,....
            'fontweight','bold','Color','r'); 
    end  
    plot(alpha_field)
  hold on
end
% Display None
if ShowElements==0 && ShowNodes==0
    cla(gcf) ;
    patch(X,Y,'w')
    plot(alpha_field)
   hold on
end

return