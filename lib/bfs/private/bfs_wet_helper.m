% Push r into queue
qt = qt + 1; % advance tail of the queue
qu(qt) = r;
%dscv(r) = true; % set root node as discovered by the BFS
wet(r) = false; % set root node as discovered by the BFS

while qt > qh
    qh = qh + 1; % advance head of the queue
    m = qu(qh); % me node; pop from head of queue
    for d = 1 : D
        n = A(m,d); % neighbour node
        if n  % check that n is not 0, i.e. a non-periodic boundary
            if wet(n)
                % n is on the surface, and undiscovered
                qt = qt + 1;  % Add n to queue
                qu(qt) = n;
                wet(n) = false; % mark n as discovered
            elseif trywet(n)
                k = K(n);
                % n is off the surface but in the ocean. Check for neutral connection
                nX = (n-1) * Xmat + 1; % = n if Xmat, 1 if not Xmat
                [s(n), t(n), x(n), success] = ntp_bottle_to_cast_mex(S(1:k,n), T(1:k,n), X(1:k,nX), s(m), t(m), x(m), X_TOL);
                if success
                    qt = qt + 1;     % Add n to queue
                    qu(qt) = n;
                    trywet(n) = false; % do not try wetting here again
                    newly_wet = newly_wet + 1;
                end
            end
        end
    end
end
