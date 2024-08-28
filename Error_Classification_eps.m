function [E_mat,E_Cumm_mat,e_per,CC,Eb_eps]=Error_Classification_eps(sig,binsize,perceptionT,sampf)

clear e111 e112 e11 e12 e13 e23 e1 e2 e3 e41 e42 e4 e_per sig.in sig.out
time = sig.time;
p_limit = perceptionT;
sign_in = sign(sig.in);

k = sig.out./max(abs(sig.in),1e-10);
k = k.*sign_in;
N = length(sig.in);

% error classification
% e1, e2, e3, e4 - logical binary value for SE, MCE, FCE, FDCE
clear i
for i=1:length(sig.in)
        if k(i)>=0
            if ((sig.in(i)>=0) && (sig.out(i)>=0))
                e111(i) = (sig.in(i)>p_limit);
                e112(i) =  (sig.out(i)>p_limit);
                e12(i) = (k(i)>1.05);
                e13(i) = (k(i)<0.95);
                e2(i) = (sig.out(i) < p_limit) && (sig.in(i) > p_limit);
                e3(i) = (sig.out(i) > p_limit) && (sig.in(i) < p_limit);
                e4(i)=0;
                e11(i) = e111(i) && e112(i);
                e23(i) = e12(i) | e13(i);
                e1(i) = e11(i) && e23(i);
            end
            if ((sig.in(i)<=0) && (sig.out(i)<=0))
                e111(i) = (sig.in(i)<-p_limit);
                e112(i) =  (sig.out(i)<-p_limit);
                e12(i) = (k(i)>1.05);
                e13(i) = (k(i)<0.95);
                e2(i) = (sig.out(i) > -p_limit) && (sig.in(i) < -p_limit);
                e3(i) = (sig.out(i) < -p_limit) && (sig.in(i) > -p_limit);
                e4(i) = 0;
                e11(i) = e111(i) && e112(i);
                e23(i) = e12(i) | e13(i);
                e1(i) = e11(i) && e23(i);
            end
        elseif  k(i)<=0
            if ((sig.in(i)>=0) && (sig.out(i)<=0))
                e2(i) = (sig.out(i) > -p_limit) && (sig.in(i) > p_limit);
                e3(i) = (sig.out(i) < -p_limit) && (sig.in(i) < p_limit);
                e1(i) = 0;
            end
            if ((sig.in(i)<=0) && (sig.out(i)>=0))
                e2(i) = (sig.out(i) < p_limit) && (sig.in(i) < -p_limit);
                e3(i) = (sig.out(i) > p_limit) && (sig.in(i) > -p_limit);
                e1(i) = 0;
            end
            e41(i) = ((sig.out(i) > p_limit) && (sig.in(i) < -p_limit));
            e42(i) = ((sig.out(i) < -p_limit) && (sig.in(i) > p_limit));
            e4(i) = e41(i) | e42(i);

        end
end
        
        
    e1 = double(e1); e2 = double(e2); e3=double(e3); e4=double(e4);
    e_per = trapz([e1' e2' e3' e4'])/N*100;

    clear E1 E2 E3 E4 E1c E2c E3c E4c
    E1 = (e1'.*abs(sig.in-sig.out));
    E2 = (e2'.*abs(sig.in-sig.out));
    E3 = (e3'.*abs(sig.in-sig.out));
    E4 = (e4'.*abs(sig.in-sig.out));
    
    E_mat = [E1 E2 E3 E4];  %instanteneous error matrix

    E1c = sum(abs(e1'.*abs(sig.in-sig.out)));
    E2c = sum(abs(e2'.*abs(sig.in-sig.out)));
    E3c = sum(abs(e3'.*abs(sig.in-sig.out)));
    E4c = sum(abs(e4'.*abs(sig.in-sig.out)));

    E_Cumm_mat = [E1c E2c E3c E4c];  % cummulative error matrix


% moving co-relation coefficient calculation

[dis_time, edge_time] = discretize(time,binsize);

for t =1:dis_time(end)
    clear find_idx
    find_idx = find(dis_time==t);
    corr_seg(t) = corr2(sig.in(find_idx),sig.out(find_idx));
    CC(find_idx(1):find_idx(end))=corr_seg(t);    
end
% moving error for different binsize normalized by bin size in seconds
E_mat(:,5)=E_mat(:,1)+1.5*E_mat(:,2)+3*E_mat(:,3)+6*E_mat(:,4);
Eb_eps(:,1) = movsum(E_mat(:,5),[sampf*0.5 0]);
Eb_eps(:,2) = movsum(E_mat(:,5),[sampf 0]);
Eb_eps(:,3) = movsum(E_mat(:,5),[sampf*2 0]);
Eb_eps(:,4) = movsum(E_mat(:,5),[sampf*5 0]);

% This is added after adding to git track
end