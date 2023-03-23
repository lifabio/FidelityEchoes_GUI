% Code inserted in the code view part of the .mlapp file

classdef Hex1_interactive_v0 < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        Label_8                       matlab.ui.control.Label
        Lamp                          matlab.ui.control.Lamp
        InitialConditons0ButtonGroup  matlab.ui.container.ButtonGroup
        Label_9                       matlab.ui.control.Label
        n0EditField                   matlab.ui.control.NumericEditField
        n0EditFieldLabel              matlab.ui.control.Label
        EditSigma                     matlab.ui.control.NumericEditField
        Label_7                       matlab.ui.control.Label
        EditRho                       matlab.ui.control.NumericEditField
        Label_6                       matlab.ui.control.Label
        n0Button                      matlab.ui.control.RadioButton
        Button                        matlab.ui.control.RadioButton
        CoherStateButton              matlab.ui.control.RadioButton
        NNButton                      matlab.ui.control.RadioButton
        Label_5                       matlab.ui.control.Label
        Label_4                       matlab.ui.control.Label
        HPLabel                       matlab.ui.control.Label
        HLabel                        matlab.ui.control.Label
        ZoomCheckBox                  matlab.ui.control.CheckBox
        NLabel                        matlab.ui.control.Label
        EditN                         matlab.ui.control.NumericEditField
        TmaxLabel                     matlab.ui.control.Label
        EditTmax                      matlab.ui.control.NumericEditField
        Label_3                       matlab.ui.control.Label
        EditEps                       matlab.ui.control.NumericEditField
        EditTheta                     matlab.ui.control.NumericEditField
        Slider_3                      matlab.ui.control.Slider
        Slider_3Label                 matlab.ui.control.Label
        EditBeta                      matlab.ui.control.NumericEditField
        EditAlpha                     matlab.ui.control.NumericEditField
        Slider_2                      matlab.ui.control.Slider
        Label_2                       matlab.ui.control.Label
        Slider                        matlab.ui.control.Slider
        Label                         matlab.ui.control.Label
        UIAxes2                       matlab.ui.control.UIAxes
        UIAxes                        matlab.ui.control.UIAxes
    end

    
    properties (Access = public)
        alpha
        beta
        theta
        eps
        Tmax
        Npt
        N
        rho
        sigma
        n0
        %Working variables
        t
        FD
        M
    end
    
    methods (Access = public)
        
        function [M FD] = compute_FidelityANDEchoes(app)
            
            alpha=app.alpha
            beta=app.beta
            theta=app.theta
            eps=app.eps
            Tmax=app.Tmax
            Npt=app.Npt
            N=app.N
            rho=app.rho;
            sigma=app.sigma;
            n0=app.n0;
            
            t=linspace(0,Tmax,Npt);
            dt=t(2)-t(1);
            %Initial conditions
            switch get(app.InitialConditons0ButtonGroup.SelectedObject,  'Text' )
                case '|φ0>= ( |1> + |2> + ... + |N> ) /√N' %where "Button1" is the name of the first radio button
                    Psi0= [0; diag(eye(N))/sqrt(N)]; 
                case '|φ0>= |α0>  Coher. State'%'Button2'
                    alpha0=rho*exp(1i*sigma);
                    Psi0= diag(zeros(N+1));
                    for n=0:N
                        Psi0(n+1)= exp(-(abs(alpha0)^2)/2)*(alpha0^n)/sqrt(factorial(n));
                    end
                case '|φ0>= |0>'%
                    Psi0= [1; diag(zeros(N))]; 
                case '|φ0>= |n0>'%
                    Psi0= diag(zeros(N+1));
                    Psi0(n0+1)=1;
                otherwise
                    Psi0= [0; diag(eye(N))/sqrt(N)]; 
            end

            a=cos(theta);
            b=sin(theta);
            for ni=0:N
                for nj=0:N
                    Kin=0;
                    a2dag=0;
                    a2=0;
                    %Kinetic energy  (alpha+beta)*(a^dag a + 1/2)
                    if ( ni==nj )
                        Kin=(alpha + beta)*(nj + 1/2);
                        Pk=(a + b)*(nj + 1/2);
                    else
                        Kin=0;
                        Pk=0;
                    end
                    %Interaction (a^dag)^2
                    if (nj==ni-2)
                        a2dag=sqrt(nj + 1)*sqrt(nj + 2);
                    else
                        a2dag=0;
                    end
                    %Interaction a^2
                    if (nj==ni+2)   
                        a2=sqrt(nj)*sqrt(nj - 1);
                    else
                        a2=0;
                    end
                    H(ni+1,nj+1)= Kin - ((alpha -beta)/2)*(a2dag + a2);
                    P(ni+1,nj+1)= Pk - ((a - b)/2)*(a2dag + a2);
                    %H(nj,ni)=H(ni,nj); %to make it more efficient together with for loop from nj=ni:Npart        
                end
            end
            H2=H + eps*P; 
            
            [V1,D1]=eig(H);
            % Lambda1(i,:)=diag(D1);
            omega1=diag(D1);
                
            [V2,D2]=eig(H2);
            omega2=diag(D2);
            
            %Expand Psi0 in the base of H
            for k=1:N+1
                c0_1(k)=V1(:,k)'*Psi0;
            end
            
            %Expand Psi0 in the base of H2
            for k=1:N+1
                c0_2(k)=V2(:,k)'*Psi0;
            end
            
            for tt=1:Npt
                %Diagonalized evolution
                Psi1=Psi0*0;
                Psi2=Psi0*0;
                Psif=Psi0*0;
                for ii=1:N+1
                    Psi1= Psi1 + exp(-1i*omega1(ii)*t(tt))*c0_1(ii)*V1(:,ii);
                    Psi2= Psi2 + exp(-1i*omega2(ii)*t(tt))*c0_2(ii)*V2(:,ii);
                end
                %Expansion of Forward propagation in H2 basis
                for k=1:N+1
                    c_fwdH2(k,1)=V2(:,k)'*Psi1;
                end
                %Backward propagation
                for ii=1:N+1
                    Psif=Psif + exp(1i*omega2(ii)*t(tt))*c_fwdH2(ii)*V2(:,ii);
                end
                
                %Fidelity and Lochsmith echoes (computed in the basis of H)
                FD(tt)=abs(Psi2'*Psi1)^2/( (Psi2'*Psi2)*(Psi1'*Psi1)  );
                M(tt)=abs(Psi0'*Psif)^2/( (Psi0'*Psi0)*(Psif'*Psif)  );
            end
        end
        
       
        
        function plot_FDandM(app)
            col=jet(8);
            col(4:5,:)=col(4:5,:)*0.9;
            vv=[2 6 4];
            cl2=col(vv,:);
            BluD=cl2(1,:);
            OrngD=cl2(2,:);
            Green=cl2(3,:);
            Violet=[113, 0, 179]/255; 
            col3=[BluD; OrngD; Green; Violet; col(8,:);col(8,:);col(7,:);col(6,:);col(5,:);col(4,:);col(3,:);col(2,:);col(1,:)];
            
            t=app.t;
            FD=app.FD;
            M=app.M;
            Tmax=app.Tmax;
            figure(app.UIFigure)
            hold(app.UIAxes, 'off')
            plot(app.UIAxes,t,FD(:),'-','Color',col3(1,:),'linewidth',1.5)
            hold(app.UIAxes, 'on')
            plot(app.UIAxes,t,M(:),'--','Color',col3(2,:),'linewidth',1.5)
            lg=legend(app.UIAxes,'F(t)','M(t)');
            xlim(app.UIAxes,[0,Tmax]);
            xlabel(app.UIAxes,'t')
            set(app.UIAxes,'Fontsize',15);
            set(lg,'box','off','Fontsize',15)
        end
    
        
        
        function plot_PhDiag(app)
            col=jet(8);
            col(4:5,:)=col(4:5,:)*0.9;
            vv=[2 6 4];
            cl2=col(vv,:);
            BluD=cl2(1,:);
            OrngD=cl2(2,:);
            Green=cl2(3,:);
            Violet=[113, 0, 179]/255; 
            col3=[BluD; OrngD; Green; Violet; col(8,:);col(8,:);col(7,:);col(6,:);col(5,:);col(4,:);col(3,:);col(2,:);col(1,:)];
            
            alpha=app.alpha
            beta=app.beta
            theta=app.theta
            eps=app.eps
            v0=linspace(-1.5,1.5,1000);
            figure(app.UIFigure)
            hold(app.UIAxes2, 'off')
            plot(app.UIAxes2,v0,0*v0,'k.','linewidth',2)
            hold(app.UIAxes2, 'on')
            plot(app.UIAxes2,0*v0,v0,'k.','linewidth',2)
            hold(app.UIAxes2, 'on')
            plot(app.UIAxes2,0,0,'.','Color',col3(2,:),'MarkerFaceColor',col3(2,:),'MarkerSize',10)
            hold(app.UIAxes2, 'on')
            plot(app.UIAxes2,alpha +eps*cos(theta),beta +eps*sin(theta),'s','Color',col3(3,:),'MarkerFaceColor',col3(3,:),'MarkerSize',7)
            hold(app.UIAxes2, 'on')
            plot(app.UIAxes2,alpha,beta,'.','Color',col3(1,:),'MarkerFaceColor',col3(1,:),'MarkerSize',20)
%             lg=legend(app.UIAxes2,'H','H + \epsilon P');
            if (app.ZoomCheckBox.Value==1)
                xlim(app.UIAxes2,[alpha-5*eps,alpha+5*eps]);
                ylim(app.UIAxes2,[beta-5*eps,beta+5*eps]);
            else
                xlim(app.UIAxes2,[-1.5,1.5]);
                ylim(app.UIAxes2,[-1.5,1.5]);
            end
            axis(app.UIAxes2,'square');
            %xlabel(app.UIAxes,'t')
            set(app.UIAxes2,'Fontsize',15);
%             set(lg,'box','off')
        end
        
        function update(app)
            app.Lamp.Color=[1, 0, 0];
            drawnow;
            %Compute Fidelity and M 
            [M FD] = compute_FidelityANDEchoes(app);
            app.FD=FD;
            app.M=M;
            %Plot
            plot_FDandM(app);
            plot_PhDiag(app);      
            app.Lamp.Color=[0, 1, 0];
            drawnow;
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            %Set initial values
            alpha=1; 
            beta=1;
            theta=pi/3;
            eps=1e-2;
            Tmax=1000;
            Npt=10000;
            N=10;
            rho=1;
            sigma=pi;
            n0=1;
            
            app.alpha=alpha;
            app.beta=beta;
            app.theta=theta;
            app.eps=eps;
            app.Tmax=Tmax;
            app.Npt=Npt;
            app.N=N;
            app.rho=rho;
            app.sigma=sigma;
            app.n0=n0;
            
            t=linspace(0,Tmax,Npt);
            dt=t(2)-t(1);
            Psi0= diag(eye(N))/sqrt(N);
            
            %Initialize bars and text
            app.Slider.Value=alpha;
            app.EditAlpha.Value=alpha;
            app.Slider_2.Value=beta;
            app.EditBeta.Value=beta;
            app.Slider_3.Value=theta;
            app.EditTheta.Value=theta/pi;
            app.EditEps.Value=eps;
            app.EditTmax.Value=Tmax;
            app.EditN.Value=N;
            app.EditRho.Value=rho;
            app.EditSigma.Value=sigma/pi;
            app.n0EditField.Value=n0;
            
            %Compute Fidelity and M with initial values
            app.Lamp.Color=[1, 0, 0];
            drawnow;
            [M FD] = compute_FidelityANDEchoes(app);
            app.t=t;
            app.FD=FD;
            app.M=M;
            %Plot
            plot_FDandM(app);
            plot_PhDiag(app);
            app.Lamp.Color=[0, 1, 0];
            drawnow;
        end

        % Value changed function: ZoomCheckBox
        function ZoomCheckBoxValueChanged(app, event)
            value = app.ZoomCheckBox.Value;
            plot_PhDiag(app);
        end

        % Value changed function: EditN
        function EditNValueChanged(app, event)
            value = app.EditN.Value;
            app.N=value;
            %Update computation and plots
            update(app);
        end

        % Value changed function: EditTmax
        function EditTmaxValueChanged(app, event)
            value = app.EditTmax.Value;
            app.Tmax=value;
            t=linspace(0,app.Tmax,app.Npt);
            app.t=t;
            %Update computation and plots
            update(app);
        end

        % Value changed function: EditEps
        function EditEpsValueChanged(app, event)
            value = app.EditEps.Value;
            app.eps=value;
            %Update computation and plots
            update(app);
        end

        % Value changing function: Slider_3
        function Slider_3ValueChanging(app, event)
            changingValue = event.Value;
            app.theta=changingValue;
            app.EditTheta.Value=changingValue/pi;
            %Update computation and plots
            update(app);
        end

        % Value changed function: EditTheta
        function EditThetaValueChanged(app, event)
            value = app.EditTheta.Value;
            app.theta=value*pi;
            app.Slider_3.Value=value*pi;
            %Update computation and plots
            update(app);
        end

        % Value changing function: Slider_2
        function Slider_2ValueChanging(app, event)
            changingValue = event.Value;
            app.beta=changingValue;
            app.EditBeta.Value=changingValue;
            %Update computation and plots
            update(app);
        end

        % Value changed function: EditBeta
        function EditBetaValueChanged(app, event)
            value = app.EditBeta.Value;
            app.beta=value;
            app.Slider_2.Value=value;
            %Update computation and plots
            update(app);
        end

        % Value changing function: Slider
        function SliderValueChanging(app, event)
            changingValue = event.Value;
            app.alpha=changingValue;
            app.EditAlpha.Value=changingValue;
            %Update computation and plots
            update(app);
        end

        % Value changed function: EditAlpha
        function EditAlphaValueChanged(app, event)
            value = app.EditAlpha.Value;
            app.alpha=value;
            app.Slider.Value=value;
            %Update computation and plots
            update(app);
        end

        % Callback function
        function InitialConditons0ButtonGroupSelectionChanged(app, event)
            selectedButton = app.InitialConditons0ButtonGroup.SelectedObject;
            %Update computation and plots
            update(app);
        end

        % Selection changed function: InitialConditons0ButtonGroup
        function InitialConditons0ButtonGroupSelectionChanged2(app, event)
            selectedButton = app.InitialConditons0ButtonGroup.SelectedObject;
            radbutton=selectedButton
            %Update computation and plots
            update(app);
        end

        % Value changed function: EditRho
        function EditRhoValueChanged(app, event)
            value = app.EditRho.Value;
            app.rho=value;
            %Update computation and plots
            update(app);
        end

        % Value changed function: EditSigma
        function EditSigmaValueChanged(app, event)
            value = app.EditSigma.Value;
            app.sigma=value*pi;
            %Update computation and plots
            update(app);
        end

        % Value changed function: n0EditField
        function n0EditFieldValueChanged(app, event)
            value = app.n0EditField.Value;
            app.n0=value;
            %Update computation and plots
            update(app);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 50 784 829];
            app.UIFigure.Name = 'UI Figure';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            xlabel(app.UIAxes, 't')
            app.UIAxes.XTickLabelRotation = 0;
            app.UIAxes.YTickLabelRotation = 0;
            app.UIAxes.ZTickLabelRotation = 0;
            app.UIAxes.Position = [59 374 683 444];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.UIFigure);
            xlabel(app.UIAxes2, '\alpha')
            ylabel(app.UIAxes2, '\beta')
            app.UIAxes2.XTickLabelRotation = 0;
            app.UIAxes2.YTickLabelRotation = 0;
            app.UIAxes2.ZTickLabelRotation = 0;
            app.UIAxes2.FontSize = 15;
            app.UIAxes2.Position = [513 145 229 230];

            % Create Label
            app.Label = uilabel(app.UIFigure);
            app.Label.HorizontalAlignment = 'right';
            app.Label.FontSize = 15;
            app.Label.Position = [47 332 25 22];
            app.Label.Text = 'α';

            % Create Slider
            app.Slider = uislider(app.UIFigure);
            app.Slider.Limits = [-1.5 1.5];
            app.Slider.ValueChangingFcn = createCallbackFcn(app, @SliderValueChanging, true);
            app.Slider.FontSize = 15;
            app.Slider.Position = [93 341 318 3];

            % Create Label_2
            app.Label_2 = uilabel(app.UIFigure);
            app.Label_2.HorizontalAlignment = 'right';
            app.Label_2.FontSize = 15;
            app.Label_2.Position = [45 263 25 22];
            app.Label_2.Text = 'ß';

            % Create Slider_2
            app.Slider_2 = uislider(app.UIFigure);
            app.Slider_2.Limits = [-1.5 1.5];
            app.Slider_2.ValueChangingFcn = createCallbackFcn(app, @Slider_2ValueChanging, true);
            app.Slider_2.FontSize = 15;
            app.Slider_2.Position = [91 272 320 3];

            % Create EditAlpha
            app.EditAlpha = uieditfield(app.UIFigure, 'numeric');
            app.EditAlpha.ValueChangedFcn = createCallbackFcn(app, @EditAlphaValueChanged, true);
            app.EditAlpha.HorizontalAlignment = 'left';
            app.EditAlpha.FontSize = 15;
            app.EditAlpha.Position = [429 332 51 22];

            % Create EditBeta
            app.EditBeta = uieditfield(app.UIFigure, 'numeric');
            app.EditBeta.ValueChangedFcn = createCallbackFcn(app, @EditBetaValueChanged, true);
            app.EditBeta.HorizontalAlignment = 'left';
            app.EditBeta.FontSize = 15;
            app.EditBeta.Position = [429 263 51 22];

            % Create Slider_3Label
            app.Slider_3Label = uilabel(app.UIFigure);
            app.Slider_3Label.HorizontalAlignment = 'right';
            app.Slider_3Label.FontSize = 15;
            app.Slider_3Label.Position = [44 192 25 22];
            app.Slider_3Label.Text = 'Θ';

            % Create Slider_3
            app.Slider_3 = uislider(app.UIFigure);
            app.Slider_3.Limits = [0 6.28318530717959];
            app.Slider_3.MajorTicks = [0 1.5707963267949 3.14159265358979 4.71238898038469 6.28318530717959];
            app.Slider_3.MajorTickLabels = {'0', 'π/2', 'π', '3π/2', '2π'};
            app.Slider_3.ValueChangingFcn = createCallbackFcn(app, @Slider_3ValueChanging, true);
            app.Slider_3.MinorTicks = [0.392699081698724 0.785398163397448 1.17809724509617 1.96349540849362 2.35619449019234 2.74889357189107 3.53429173528852 3.92699081698724 4.31968989868597 5.10508806208341 5.49778714378214 5.89048622548086];
            app.Slider_3.FontSize = 15;
            app.Slider_3.Position = [90 201 320 3];

            % Create EditTheta
            app.EditTheta = uieditfield(app.UIFigure, 'numeric');
            app.EditTheta.ValueChangedFcn = createCallbackFcn(app, @EditThetaValueChanged, true);
            app.EditTheta.HorizontalAlignment = 'left';
            app.EditTheta.FontSize = 15;
            app.EditTheta.Position = [429 192 51 22];

            % Create EditEps
            app.EditEps = uieditfield(app.UIFigure, 'numeric');
            app.EditEps.ValueChangedFcn = createCallbackFcn(app, @EditEpsValueChanged, true);
            app.EditEps.HorizontalAlignment = 'left';
            app.EditEps.FontSize = 15;
            app.EditEps.Position = [80 124 51 22];

            % Create Label_3
            app.Label_3 = uilabel(app.UIFigure);
            app.Label_3.FontSize = 15;
            app.Label_3.Position = [59 124 25 22];
            app.Label_3.Text = 'ε';

            % Create EditTmax
            app.EditTmax = uieditfield(app.UIFigure, 'numeric');
            app.EditTmax.ValueChangedFcn = createCallbackFcn(app, @EditTmaxValueChanged, true);
            app.EditTmax.HorizontalAlignment = 'left';
            app.EditTmax.FontSize = 15;
            app.EditTmax.Position = [242 124 51 22];

            % Create TmaxLabel
            app.TmaxLabel = uilabel(app.UIFigure);
            app.TmaxLabel.FontSize = 15;
            app.TmaxLabel.Position = [194 124 43 22];
            app.TmaxLabel.Text = 'Tmax';

            % Create EditN
            app.EditN = uieditfield(app.UIFigure, 'numeric');
            app.EditN.ValueChangedFcn = createCallbackFcn(app, @EditNValueChanged, true);
            app.EditN.HorizontalAlignment = 'left';
            app.EditN.FontSize = 15;
            app.EditN.Position = [429 124 51 22];

            % Create NLabel
            app.NLabel = uilabel(app.UIFigure);
            app.NLabel.HorizontalAlignment = 'center';
            app.NLabel.FontSize = 15;
            app.NLabel.Position = [398 124 25 22];
            app.NLabel.Text = 'N';

            % Create ZoomCheckBox
            app.ZoomCheckBox = uicheckbox(app.UIFigure);
            app.ZoomCheckBox.ValueChangedFcn = createCallbackFcn(app, @ZoomCheckBoxValueChanged, true);
            app.ZoomCheckBox.Text = 'Zoom';
            app.ZoomCheckBox.FontSize = 15;
            app.ZoomCheckBox.Position = [569 124 60 22];

            % Create HLabel
            app.HLabel = uilabel(app.UIFigure);
            app.HLabel.FontSize = 15;
            app.HLabel.Position = [700 132 25 22];
            app.HLabel.Text = 'H';

            % Create HPLabel
            app.HPLabel = uilabel(app.UIFigure);
            app.HPLabel.FontSize = 15;
            app.HPLabel.Position = [700 111 50 22];
            app.HPLabel.Text = 'H + εP';

            % Create Label_4
            app.Label_4 = uilabel(app.UIFigure);
            app.Label_4.HorizontalAlignment = 'center';
            app.Label_4.VerticalAlignment = 'top';
            app.Label_4.FontSize = 40;
            app.Label_4.FontColor = [0 0.502 1];
            app.Label_4.Position = [667 119 25 50];
            app.Label_4.Text = '•';

            % Create Label_5
            app.Label_5 = uilabel(app.UIFigure);
            app.Label_5.HorizontalAlignment = 'center';
            app.Label_5.VerticalAlignment = 'top';
            app.Label_5.FontSize = 20;
            app.Label_5.FontColor = [0.451 0.902 0.451];
            app.Label_5.Position = [658 91 43 46];
            app.Label_5.Text = '■';

            % Create InitialConditons0ButtonGroup
            app.InitialConditons0ButtonGroup = uibuttongroup(app.UIFigure);
            app.InitialConditons0ButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @InitialConditons0ButtonGroupSelectionChanged2, true);
            app.InitialConditons0ButtonGroup.Title = 'Initial Conditons |φ0>:';
            app.InitialConditons0ButtonGroup.FontSize = 15;
            app.InitialConditons0ButtonGroup.Position = [31 13 719 91];

            % Create NNButton
            app.NNButton = uiradiobutton(app.InitialConditons0ButtonGroup);
            app.NNButton.Text = '|φ0>= ( |1> + |2> + ... + |N> ) /√N';
            app.NNButton.FontSize = 13;
            app.NNButton.Position = [10 36 255 27];
            app.NNButton.Value = true;

            % Create CoherStateButton
            app.CoherStateButton = uiradiobutton(app.InitialConditons0ButtonGroup);
            app.CoherStateButton.Text = '|φ0>= |α0>  Coher. State';
            app.CoherStateButton.FontSize = 13;
            app.CoherStateButton.Position = [264 38 165 22];

            % Create Button
            app.Button = uiradiobutton(app.InitialConditons0ButtonGroup);
            app.Button.Text = '|φ0>= |0>';
            app.Button.FontSize = 13;
            app.Button.Position = [472 38 78 22];

            % Create n0Button
            app.n0Button = uiradiobutton(app.InitialConditons0ButtonGroup);
            app.n0Button.Text = '|φ0>= |n0>';
            app.n0Button.FontSize = 13;
            app.n0Button.Position = [588 38 85 22];

            % Create Label_6
            app.Label_6 = uilabel(app.InitialConditons0ButtonGroup);
            app.Label_6.HorizontalAlignment = 'right';
            app.Label_6.FontSize = 13;
            app.Label_6.Position = [258 9 27 22];
            app.Label_6.Text = '|α0|';

            % Create EditRho
            app.EditRho = uieditfield(app.InitialConditons0ButtonGroup, 'numeric');
            app.EditRho.ValueChangedFcn = createCallbackFcn(app, @EditRhoValueChanged, true);
            app.EditRho.FontSize = 15;
            app.EditRho.Position = [292 9 37 22];

            % Create Label_7
            app.Label_7 = uilabel(app.InitialConditons0ButtonGroup);
            app.Label_7.HorizontalAlignment = 'right';
            app.Label_7.FontSize = 13;
            app.Label_7.Position = [352 9 25 22];
            app.Label_7.Text = 'σ';

            % Create EditSigma
            app.EditSigma = uieditfield(app.InitialConditons0ButtonGroup, 'numeric');
            app.EditSigma.ValueChangedFcn = createCallbackFcn(app, @EditSigmaValueChanged, true);
            app.EditSigma.FontSize = 15;
            app.EditSigma.Position = [383 9 37 22];

            % Create n0EditFieldLabel
            app.n0EditFieldLabel = uilabel(app.InitialConditons0ButtonGroup);
            app.n0EditFieldLabel.HorizontalAlignment = 'right';
            app.n0EditFieldLabel.FontSize = 13;
            app.n0EditFieldLabel.Position = [597 9 25 22];
            app.n0EditFieldLabel.Text = 'n0';

            % Create n0EditField
            app.n0EditField = uieditfield(app.InitialConditons0ButtonGroup, 'numeric');
            app.n0EditField.ValueChangedFcn = createCallbackFcn(app, @n0EditFieldValueChanged, true);
            app.n0EditField.FontSize = 15;
            app.n0EditField.Position = [628 9 37 22];

            % Create Label_9
            app.Label_9 = uilabel(app.InitialConditons0ButtonGroup);
            app.Label_9.HorizontalAlignment = 'center';
            app.Label_9.FontSize = 13;
            app.Label_9.Position = [419 9 25 22];
            app.Label_9.Text = 'π';

            % Create Lamp
            app.Lamp = uilamp(app.UIFigure);
            app.Lamp.Position = [740 785 20 20];

            % Create Label_8
            app.Label_8 = uilabel(app.UIFigure);
            app.Label_8.HorizontalAlignment = 'center';
            app.Label_8.FontSize = 15;
            app.Label_8.Position = [480 192 25 22];
            app.Label_8.Text = 'π';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Hex1_interactive_v0

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
