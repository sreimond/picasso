

function normals2regtools(file_path,reg_type)
[N,n,x,n_obs,n_para,lPl] = read_groops_normals( file_path );
out_x = [file_path '-x_regtools.txt'];
out_sigmax = [file_path '-sigmax_regtools.txt'];
A = N;
b = n;
[U,s,V] = csvd(A);
if strcmpi(reg_type,'lcurve')
    lambda = l_curve(U,s,b);
elseif strcmpi(reg_type,'gcv')
    lambda = gcv(U,s,b);
end
R = eye(size(N))*lambda;
x = (N+R)\n;
rss = x'*N*x + lPl - 2*n'*x;
sig0 = rss/(n_obs-n_para);
QQ = inv(N+R);
Q = QQ*N*QQ;
Q = Q * sig0;
%dlmwrite(sig0,'/home/sreimond/tbd_sig0.txt')
%dlmwrite([n_obs,n_para,n_obs-n_para],'/home/sreimond/tbd_n.txt')

%fprintf('N(1,1) = %+.5e\n',N(1,1))
%fprintf('R(1,1) = %+.5e\n',R(1,1))
%fprintf('Any N negative?: %d\n',any(N(:)<0))
%fprintf('Any R negative?: %d\n',any(R(:)<0))
%fprintf('Any (N+R) negative?: %d\n',any((N(:)+R(:))<0))
%QQ1 = inv(N+R);
%QQ = QQ1*N*QQ1;
%QQ = QQ1*QQ1;
%QQ = inv(N+R);
%fprintf('Any isinf(inv(N + R))?: %d\n',any(isinf(QQ(:))))
%fprintf('Any diag(inv(N+R)) negative?: %d\n',any(diag(QQ)<0))
%fprintf('Q(1,1) = %+.5e\n',Q(1,1))
%fprintf('sig0 = %.5f\n',sig0)
%fprintf('n_obs = %.5f\n',n_obs)
%fprintf('n_para = %.5f\n',n_para)
sigma_x = abs(diag(Q)).^0.5;
%sigma_x = diag(Q).^0.5;
output_ascii = zeros(numel(x),1);
output_ascii(:) = x;
dlmwrite(out_x,output_ascii,'delimiter',' ','precision','%.16e');
output_ascii(:) = sigma_x;
dlmwrite(out_sigmax,output_ascii,'delimiter',' ','precision','%.16e');
end

