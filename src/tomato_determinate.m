function [Wnow,WFnow,LAInow,Nnow,Wprev,WFprev,LAIprev, Nprev,Ns,Nl,Nf,Rc_prev,POL,POF] = ...
    tomato_determinate(plantingdensity_tomato, Nm, Nprev,LAIprev, LAImax,...
    E,D,LFmax,K,m,Qe,WFprev,Wprev,Q10,RL,RF,Tcrit,NFF,...
    DLI_inside,th,sunrise,sunset,day,CO2,Ns,Nl,Nf,day1,Thours,FPN,nL,SLAmin,SLAmax,betaT,betaC,...
    FRPET,FRSTM,TPL,tom_a1,tom_a2,tom_b1,tom_b2,tom_c1,tom_c2,rFmax,rLmax)


PPFD=DLI_inside*1e6/((sunset-sunrise)*3600);

    T=th;
    th_temp=Thours;
    Tdaytime=mean(th_temp(ceil(sunrise):ceil(sunset),1));


T1=T;
Tdaytime1=Tdaytime(1,1);


SLA=(SLAmin+(SLAmax-SLAmin)*exp(-0.471*DLI_inside))/((1+(betaT*(24-T1)))*(1+(betaC*(CO2-350))));
if(SLA<SLAmin)
    SLA=SLAmin;
end
if(SLA>SLAmax)
    SLA=SLAmax;
end

FCO2_1=1+(0.0003*(CO2-350));

POL=tom_a1*tom_b1*exp(-tom_b1*(day-tom_c1))*exp(-exp(-tom_b1*(day-tom_c1)));
POF=tom_a2*tom_b2*exp(-tom_b2*(day1-tom_c2))*exp(-exp(-tom_b2*(day1-tom_c2)));
fnT2=1+(0.0281*(Tdaytime1-28));
POF=POF/fnT2;
POL=POL/fnT2;



%Jones tomgro 1991
if(T1<12)
    rF=0;
    rL=0;
else
    if(T1>12&&T1<=Tcrit)
        rL=(0.0005*T1) - 0.0014;
        rF= ((0.0017*T1) - 0.0147);
        if(rF>rFmax)
            rF=rFmax;
        end

        if(rL>rLmax)
            rL=rLmax;
        end
    else
        rL = -0.0009*T1 + 0.0433;
        rF = ((-0.0021*T1) + 0.1067);
        if(rL<0)
            rL=0;
        end
        if(rF<0)
            rF=0;
        
        end
    end
end

if(T1<=12)
    fnT1=0;
else
    if(T1>12&&T1<=28)
        fnT1=1+(0.0281*(T1-28));
    else
        fnT1=1-(0.0455*(T1-28));
        if(fnT1<0)
            fnT1=0;
        end
    end
end



if(day==1)
    Ns=zeros(nL,1);
    Nl=zeros(nL,1);
    Nf=zeros(nL,1);
end
for i=1:nL
    if(i==1)
    Ns(i,1)=Ns(i,1)+(Nm*fnT1*FCO2_1*plantingdensity_tomato)-...
    (rL*FCO2_1*nL*Ns(i,1));
    Nl(i,1)=Nl(i,1)+(Nm*fnT1*FCO2_1*plantingdensity_tomato/(1+TPL))-...
    (rL*FCO2_1*nL*Nl(i,1));
    else
        if(i<=nL-1)
         Ns(i,1)=Ns(i,1)+(rL*FCO2_1*nL*Ns(i-1,1))-...
             (rL*FCO2_1*nL*Ns(i,1));
         Nl(i,1)=Nl(i,1)+(rL*FCO2_1*nL*Nl(i-1,1))-...
             (rL*FCO2_1*nL*Nl(i,1));
        else
         Ns(i,1)=Ns(i,1)+(rL*FCO2_1*nL*Ns(i-1,1));
         Nl(i,1)=Nl(i,1)+(rL*FCO2_1*nL*Nl(i-1,1));
        end

    end
    if(i==1)
        if(Nprev>=NFF)
            Nf(i,1)=Nf(i,1)+(Nm*fnT1*FCO2_1*plantingdensity_tomato*FPN)-(rF*FCO2_1*nL*Nf(i,1));
        else
            Nf(i,1)=0;
            POF=0;
        end
    else
        if(Nprev>=NFF&&i<=nL-1)
            Nf(i,1)=Nf(i,1)-(rF*FCO2_1*nL*Nf(i,1))+(rF*FCO2_1*nL*Nf(i-1,1));
        end
        if(Nprev>=NFF&&i==nL)
            Nf(i,1)=Nf(i,1)+(rF*FCO2_1*nL*Nf(i-1,1));
        end
        if(Nprev<NFF)
            Nf(i,1)=0;
            POF=0;
        end
    end
    
end


POL_matrix=POL*fnT1*FCO2_1;
POF_matrix=POF*fnT1*FCO2_1;
dAlpdt=Nl.*POL_matrix;
Ldem_matrix=(1+FRPET)*(dAlpdt/SLA);
Ldem=sum(Ldem_matrix); %g dw/m2-day
LAI_dem=LAIprev+(Ldem*SLA);
Sdem_matrix=Ldem_matrix.*FRSTM.*(Ns./Nl);
Sdem=sum(Sdem_matrix);
Fdem_matrix=POF_matrix.*Nf;
Fdem=sum(Fdem_matrix);
if(Fdem>0)
    stophere=1;
end
Demand=Ldem+Sdem+Fdem;
%% Temperature and node dependent properties
% Truss (flowering) and hence nodes (assumed)have
%appearance rate that is equally related to day and
%night temperature and therefore responds to the
%24 hour mean temperature. No effects of fruit load
%leaf removal or plant density on truss appearance
%rate have been observed (Tomatoes book page 71)
%Sink-source ratio has no significant influence on
%truss appearance rate. Air CO2 concentration and
%air humidity hardly affect truss appearance rate

temp = [0,9,12,28,50];
%This would be daytime temperature
temp2 = [0,9,12,15,21,28,35,50];
nodenumber = [0,1.7,4.3,5.9,8,10,12,15,20,100];
fRN = [0.15,0.15,0.16,0.16,0.18,0.18,0.19,0.2,0.2,0.2];


%Adaptation of the CROPGRO model ro simulate the growth
%of field grown tomato
%J.M.S. Scholberg, K.J. Boote, J.W. Jones and B.L. McNeal
temp_size=size(temp);
temp_size=temp_size(2);

temp2_size=size(temp2);
temp2_size=temp2_size(2);

nodenumber_size=size(nodenumber);
nodenumber_size=nodenumber_size(2);

    if(Tdaytime1<10||Tdaytime1>45)
        PGREDT1=0;
    else
        %Photosynthetic and respiratory characterization of field grown
        %tomatoes Jorge A. Bolanos & Theodore C. Hsiao
        PGREDT1=-0.0022*(Tdaytime1^2) + 0.1182*(Tdaytime1)-0.6401;
    end



    ctrl=0;
    if(Nprev>nodenumber(nodenumber_size))
        fRN1=fRN(nodenumber_size);
        ctrl=1;
    else
        for i=1:nodenumber_size
            if(Nprev==nodenumber(i))
                fRN1=fRN(i);
                ctrl=1;
                break;
            
            end
        end
    end
    if(ctrl==0)
        for i=1:nodenumber_size-1
            if(Nprev>nodenumber(i)&&Nprev<nodenumber(i+1))
                fRN1=(((fRN(i+1)-fRN(i))/(nodenumber(i+1)-nodenumber(i)))*(Nprev-nodenumber(i)...
                    ))+fRN(i);
            
            end
        end
    end

%% State variable  1
    dNdt = Nm*fnT1*FCO2_1;
    Nnow=dNdt+Nprev;
    if(Nnow>=NFF)
        Nnow=NFF;
    end
%% State variable 2

    %dLAIdt = plantingdensity_tomato*delta*lambdaTd1*(exp(beta*(Nnow-Nb))/(1+exp(beta*(Nnow-Nb))))*dNdt;
    %LAInow = dLAIdt + LAIprev;



%% State variable 3

    Pg = (D*LFmax*PGREDT1/K)*log((((1-m)*LFmax)+(Qe*K*PPFD(1,1)))/(((1-m)*LFmax)+(Qe*K*PPFD(1,1)*...
        exp(-K*LAIprev))));
    Rm = (Q10^(0.1*(T1-20)))*((RL*(Wprev))+(RF*WFprev));

    GRnet = E*(Pg-Rm)*(1-fRN1);

Rc=GRnet/Demand;
if(Demand==0)
    Rc=0;
end
if(Nprev<=5)
    if(Rc>1)
        Rc=1;
        LAInow=LAI_dem;
    else
        Ldem=Ldem*Rc; %g dw/m2-day
        LAInow=LAIprev+(Ldem*SLA);
        Sdem=Sdem*Rc;
        Fdem=Fdem*Rc;
%         Nl=Nl*Rc;
%         Ns=Ns*Rc;
%         Nf=Nf*Rc;
        Demand=Ldem+Sdem+Fdem;
    end
else
%         if(Rc>1)
%         Rc=1;
%         LAInow=LAI_dem;
%         else
    Ldem=Ldem*Rc; %g dw/m2-day
    LAInow=LAIprev+(Ldem*SLA);
    Sdem=Sdem*Rc;
    Fdem=Fdem*Rc;
    Demand=Ldem+Sdem+Fdem;
end
    if (LAInow >LAImax)
        LAInow = LAImax;
    end
    WFnow = WFprev + Fdem;


%% Tabulating parameters
    Wnow=Wprev+Ldem+Sdem;
    Wprev=Wnow;
    WFprev=WFnow;
    LAIprev=LAInow;
    Nprev=Nnow;
    Rc_prev=Rc;
    if(day==33)
        stophere=1;
    end

end
