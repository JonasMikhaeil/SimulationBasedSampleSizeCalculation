library(shiny)
library(broom)
library(ggplot2)
library(dplyr)
library(LaplacesDemon)
library(sandwich)
library(patchwork)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel(
    h1("Simulation based Sample Size Calculations", align = "center")
  ),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      h4("Simulation Parameters",align = "center"),
      fluidRow(
     column(6,
     checkboxInput("treat", "Turn treatment on", value = FALSE),
     sliderInput(inputId = "comp_prob",
                 label = "P(Compliance)",
                 min = 0,
                 max = 1,
                 value = 0.8,step = 0.05),),
     column(6,
            numericInput("num", 
                         h5("# of simulated units"), 
                         value = 10000),
     ),),
     fluidRow(
      column(6,
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "mean_col",
                  label = "mean Columbia",
                  min = -1,
                  max = 1,
                  value = -0.2,step=0.1,), #0.08
      sliderInput(inputId = "mean_all",
                  label = "mean Allen",
                  min = 0,
                  max = 1,
                  value = 0.2), #0.0
      sliderInput(inputId = "mean_law",
                  label = "mean Lawrence",
                  min = 0,
                  max = 1,
                  value = 0.2), #0.29
      sliderInput(inputId = "mean_cor",
                  label = "mean Cornell",
                  min = -1,
                  max = 1,
                  value = -0.7,step=0.1), #0
      sliderInput(inputId = "mean_lm",
                  label = "mean LowMan",
                  min = -1,
                  max = 1,
                  value = -0.1,step=0.1), #0.25
      sliderInput(inputId = "mean_qu",
                  label = "mean Queens",
                  min = 0,
                  max = 1,
                  value = 0.5), #0.24
      sliderInput(inputId = "sd",
                  label = "standard deviation",
                  min = 0,
                  max = 1,
                  value = 0.59), #0.59
      
      ),
    
      column(6,
      sliderInput(inputId = "alpha_col",
                  label = "ML miscalibration Columbia",
                  min = -3,
                  max = 3,
                  value = 0.1,step = 0.06),
      sliderInput(inputId = "alpha_all",
                  label = "ML miscalibration Allen",
                  min = -3,
                  max = 3,
                  value = -0.35,step = 0.05),
      sliderInput(inputId = "alpha_law",
                  label = "ML miscalibration Lawrence",
                  min = -3,
                  max = 3,
                  value = 0.5,step = 0.05),
      sliderInput(inputId = "alpha_cor",
                  label = "ML miscalibration Cornell",
                  min = -3,
                  max = 3,
                  value = -0.1,step = 0.05),
      sliderInput(inputId = "alpha_lm",
                  label = "ML miscalibration LowMan",
                  min = -3,
                  max = 3,
                  value = -0.35,step = 0.05),
      sliderInput(inputId = "alpha_qu",
                  label = "ML miscalibration Queens",
                  min = -3,
                  max = 3,
                  value = -0.15,step = 0.05),
      actionButton("power_calc", "Sample Size Calculation"),
      
      
      ),
     ),
    ),
     
     
   
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "depPlot"),
      column(7,plotOutput(outputId = "compPlot")),
      column(5,plotOutput(outputId = "distPlot")),
      tableOutput(outputId = "treatment_effect"),
      plotOutput(outputId = "plot")
      
      
    )
  )
)


# Define server logic required to draw a histogram ----
server <- function(input, output) {
  probabiltiy_shd_given_riskScore <- function(risk_score,hospital){
    ### assume score is calibrated
    hosp<-c("columbia","allen","lawrence","cornell","low man","queens")
    
    alpha_hosp <- c(input$alpha_col,input$alpha_all,input$alpha_law,input$alpha_cor,input$alpha_lm,input$alpha_qu)
    alpha <- alpha_hosp[which(hospital==hosp)]
    invlogit(alpha+logit(risk_score))
  }
  
  p_echo_givenRS_ShD_treatment <- function(risk_score , shd , treatment_indicator){
    X <- input$comp_prob
    if (treatment_indicator ==0){
      if(shd == 1){
        p<-0.4+(0.8-0.4)/0.3*(risk_score-0.6)
      } 
      else{
        p<-0.4+(0.8-0.4)/0.3*(risk_score-0.6)
      } 
    }
    else{
      if(shd == 1){
        RZ <- rbinom(1,1,X)
        #p<-1^RZ*(0.4+(0.8-0.4)/0.3*(risk_score-0.6))^(1-RZ)
        p<-X + (1-X)*(0.4+(0.8-0.4)/0.3*(risk_score-0.6))
      }
      else{
        RZ <- rbinom(1,1,X)
        p<-X + (1-X)*(0.4+(0.8-0.4)/0.3*(risk_score-0.6))
        #p<-1^RZ*(0.4+(0.8-0.4)/0.3*(risk_score-0.6))^(1-RZ)
        
      }
    }
    p
  }
  simulate_data <- function(N){
    K <-6
    proportions_assigned_to_hospital <- c(0.2327044,0.1228834,0.1482825,0.2005322,0.1064344,0.1891630)
    hosp<-c("columbia","allen","lawrence","cornell","low man","queens")
    hosp_mu <- c(input$mean_col,input$mean_all,input$mean_law,input$mean_cor,input$mean_lm,input$mean_qu) #c(0.005,0.06,0.12,0.02,0.03,0.11)
    hosp_var <- rep(input$sd^2,6)#c(0.005,0.06,0.12,0.02,0.03,0.11)
    
    N_K <- round(proportions_assigned_to_hospital*N)
    
    
    risk_score_vector <- numeric(0)
    
    # Generate risk scores for each hospital
    for (k in 1:length(N_K)) {
      mock <- c(0.7,0.7,0.5,0.5,0.7,0.7)
      risk_score_vector <- c(risk_score_vector,truncnorm::rtruncnorm(N_K[k],a=0.6,b=1.,hosp_mu[k],hosp_var[k])) #rep(mock[k],N_K[k]))#
    }
    
    # Create the data frame
    risk_score <- data.frame(
      hospital = rep(hosp, times = N_K),
      risk_score = risk_score_vector
    )
    
    data <- risk_score %>% 
      rowwise() %>% 
      mutate(p_shd = probabiltiy_shd_given_riskScore(risk_score,hospital)) %>%
      mutate(shd = rbinom(1, 1, p_shd)) %>%
      ungroup() 
    #make this restricted (cornell & columbia)
    if(input$treat==TRUE){
      treatment_indicator <- sample(1:K,round(K/2),replace = FALSE)
      while (which(hosp == "columbia") %in% treatment_indicator == which(hosp == "cornell") %in% treatment_indicator) {
        treatment_indicator <- sample(1:K, round(K/2), replace = FALSE)
      }}
    else treatment_indicator <-c()
    
    data$treatment <- ifelse(data$hospital %in% hosp[treatment_indicator], 1, 0)
    
    data <- data %>% 
      rowwise() %>% 
      mutate(p_echo = p_echo_givenRS_ShD_treatment(risk_score,shd,treatment))  %>% 
      mutate(echo=rbinom(1, 1, p_echo)) %>% 
      ungroup() 
    data$diagnosis <- ifelse(data$shd==1 & data$echo==1, 1, 0)
    
    data
  }
 
  
  data_reactive <- reactive({
    simulate_data(input$num)
  })
  output$depPlot <- renderPlot({
    data<-data_reactive()
    data$treatment <- factor(data$treatment, levels=c("0", "1"), labels=c("control", "treatment"))
    plot1<- ggplot(data,aes(x=risk_score,y=p_echo,color=treatment))+
      scale_color_manual(values=c("black", "steelblue")) +
      geom_line()+
      theme_minimal()+
      theme(legend.position = c(0.8,0.3),
            legend.title=element_blank())+
      ylab("P(ECHO | Risk Score)")
    plot2 <- 
      ggplot(data,aes(x=risk_score,y=p_shd,group=hospital))+
      geom_line()+
      theme_minimal() +
      ylab("P(SHD | Risk Score)")
    plot <- plot1 + plot2 & xlab(NULL) & theme(plot.margin = margin(5.5, 5.5, 0, 5.5))
    wrap_elements(panel = plot) +
      labs(tag = "Risk Score") +
      theme(
        plot.tag = element_text(size = rel(1)),
        plot.tag.position = "bottom"
      )
  })
  output$distPlot <- renderPlot({
    data<-data_reactive()
   
    N <- c(rep(0.65,1550),rep(0.75,1226),rep(0.85,900),rep(0.95,433)) 
   
    hist_dat <-data.frame(
      value=c(N,data$risk_score),
      cat= rep(c("h1", "h2"), times = c(length(N),length(data$risk_score))))
    ggplot(hist_dat, aes(x = value, fill = cat,y=..density..)) +
    geom_histogram(position = "identity", alpha = 0.5, breaks = c(0.6,0.7,0.8,0.9,1)) +
    scale_fill_manual( values=c("black", "steelblue")) +
    theme_minimal() +
    labs(title = "", x = "Risk Score", y = "")+
    theme(axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),legend.position = "none")
    
    #ggplot(hist_dat, aes(x=value,fill=cat,y=..density..)) + 
      #geom_histogram(position = "identity",  breaks = c(0.6,0.7,0.8,0.9,1), size=1) +
      #scale_fill_manual(values=c("black")) +
      #geom_histogram(position = "identity", breaks = c(0.6,0.7,0.8,0.9,1)) +
      #scale_fill_manual(values=c("white")) +
      #theme_minimal() +
      #labs(title = "", x = "Risk Score", y = "")+
      #theme(axis.text.y=element_blank(),
            #axis.ticks.y=element_blank(),legend.position = "none")
    
  })
  
 
  output$treatment_effect <- renderTable(
    {data<-data_reactive()
    VaryAClus <- data %>%
      group_by(treatment,hospital) %>%
      summarize(
        MeanyAClus = mean(diagnosis)
      ) %>%
      ungroup(.) %>%
      group_by(treatment) %>%
      summarize(
        VaryAClus = var(MeanyAClus)
      ) 
    # Total variance
    VaryA <- data %>%
      group_by(treatment) %>%
      summarize(
        VaryATot = var(diagnosis)
      )
      averages_by_treatment <- data %>%
        group_by(treatment) %>%
        summarise(
          diagnosis = mean(diagnosis),
          echo = mean(echo),
          "diagnosis per echo" =mean(diagnosis)/mean(echo),
          PPV = mean(shd)
        )   %>%
        left_join(VaryAClus,by='treatment') %>%
        left_join(VaryA,by='treatment') %>%
        mutate(
          ICC =  VaryAClus/VaryATot)
        
      # Print the results
      print(averages_by_treatment|> select("treatment","ICC","diagnosis","echo","diagnosis per echo","PPV"))
    }
  )
  
  
  output$compPlot <- renderPlot({
    data<-data_reactive()

    hosp<-c("columbia","allen","lawrence","cornell","low man","queens")
    real_data <- data.frame(hospital=hosp,shd_rate=c(299/(111+299),207/(207+84),301/(301+51),282/(135+282),166/(72+166),372/(372+111)))
    rate_of_detected_shd <- data |>filter(echo == 1) |> aggregate(diagnosis~hospital,mean)
    #colnames(rate_of_detected_shd)[3] <- "rate of shd in patients who \n received an echo (usual care)"
    plot1 <- ggplot(real_data,aes(x=hospital,y=shd_rate))+
      geom_point() + 
      geom_point(data=rate_of_detected_shd,aes(x=hospital,y=diagnosis),col="steelblue")+
      theme_minimal() +
      labs(title = "", x = "Hospital", y = "P(SHD | ECHO received)")+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    
    
    
    hosp<-c("columbia","allen","lawrence","cornell","low man","queens")
    real_data <- data.frame(hospital=hosp,shd_rate=c(299/(111+299),207/(207+84),301/(301+51),282/(135+282),166/(72+166),372/(372+111))
                            ,prob_echo =c(1-0.466,1-.427,1-.426,1-.497,1-.459,1-.382) )
    rate_of_echo <- data  |> aggregate(echo~hospital,mean)
    #colnames(rate_of_detected_shd)[3] <- "rate of shd in patients who \n received an echo"
    plot2 <- ggplot(real_data,aes(x=hospital,y=prob_echo))+
      geom_point() + 
      geom_point(data=rate_of_echo,aes(x=hospital,y=echo),col="steelblue")+
      theme_minimal() +
    labs(title = "", x = "Hospital", y = "P(ECHO)")+
    annotate(geom="text", x="queens", y=0.515, label="Simulated",color="steelblue")+
    annotate(geom="text", x="queens", y=0.50, label="Data",color="black")+
      theme(text = element_text(size=13))
    plot2 / plot1 & xlab(NULL) & theme(plot.margin = margin(5.5, 5.5, 0, 5.5))
    
    
  })
  
 
  powerVals <- eventReactive(input$power_calc, {
    simulate_data <- function(N){
      K <-6
      proportions_assigned_to_hospital <- c(0.2327044,0.1228834,0.1482825,0.2005322,0.1064344,0.1891630)
      hosp<-c("columbia","allen","lawrence","cornell","low man","queens")
      hosp_mu <- c(input$mean_col,input$mean_all,input$mean_law,input$mean_cor,input$mean_lm,input$mean_qu) #c(0.005,0.06,0.12,0.02,0.03,0.11)
      hosp_var <- rep(input$sd^2,6)#c(0.005,0.06,0.12,0.02,0.03,0.11)
      
      N_K <- round(proportions_assigned_to_hospital*N)
      
      
      risk_score_vector <- numeric(0)
      
      # Generate risk scores for each hospital
      for (k in 1:length(N_K)) {
        mock <- c(0.7,0.7,0.5,0.5,0.7,0.7)
        risk_score_vector <- c(risk_score_vector,truncnorm::rtruncnorm(N_K[k],a=0.6,b=.9,hosp_mu[k],hosp_var[k])) #rep(mock[k],N_K[k]))#
      }
      
      # Create the data frame
      risk_score <- data.frame(
        hospital = rep(hosp, times = N_K),
        risk_score = risk_score_vector
      )
      
      data <- risk_score %>% 
        rowwise() %>% 
        mutate(p_shd = probabiltiy_shd_given_riskScore(risk_score,hospital)) %>%
        mutate(shd = rbinom(1, 1, p_shd)) %>%
        ungroup() 
      #make this restricted (cornell & columbia)
      
      treatment_indicator <- sample(1:K,round(K/2),replace = FALSE)
      while (which(hosp == "columbia") %in% treatment_indicator == which(hosp == "cornell") %in% treatment_indicator) {
        treatment_indicator <- sample(1:K, round(K/2), replace = FALSE)
      }
     
      
      data$treatment <- ifelse(data$hospital %in% hosp[treatment_indicator], 1, 0)
      
      data <- data %>% 
        rowwise() %>% 
        mutate(p_echo = p_echo_givenRS_ShD_treatment(risk_score,shd,treatment))  %>% 
        mutate(echo=rbinom(1, 1, p_echo)) %>% 
        ungroup() 
      data$diagnosis <- ifelse(data$shd==1 & data$echo==1, 1, 0)
      
      data
    }
    calculate_power <- function(N,clustering=TRUE){
      
      alpha <- 0.05
      sig_level <- alpha/2
      
      rep <- seq(1:100)
      
      # Export necessary functions and libraries to the parallel environment
      
      power <-sapply(rep,function(r){
        data <- simulate_data(N)
        #mod <- lm(diagnosis~treatment+risk_score,data=data)
        
        #var_te <- car::hccm(mod,type="hc2")["treatment","treatment"]
        #coef <- tidy(mod)
        ###using EHW standard errors (see (Lei and Ding, 2021))
        #if(clustering==TRUE){
          #VCovClusterHC1 <- sandwich::vcovCL(mod,cluster=data$hospital,type="HC1",sandwich=TRUE)
          #print(VCovClusterHC1)
          #var_te <-VCovClusterHC1["treatment","treatment"]
          #if(abs(coef$estimate[2]) - qnorm(1-sig_level) * sqrt(var_te)  > 0){r = 1}
          #else{r=0}
        #}
        if(clustering==TRUE){
          mod <- lme4::lmer(diagnosis~1+treatment+risk_score + (1| hospital),data=data)
          
          te <- fixef(mod)["treatment"]
          te.se <- se.fixef(mod)["treatment"]
          if(abs(te) - qnorm(1-sig_level) * te.se  > 0){r = 1}
          else{r=0}
        }
        else{
          mod <- lm(diagnosis~treatment+risk_score,data=data)
          coef <- tidy(mod)
          var_te <- sandwich::vcovHC(mod, type = "HC1")["treatment","treatment"]
          sd_te<- coef$std.error[2]#sqrt(var_te)
          if(abs(coef$estimate[2]) - qnorm(1-sig_level) * sqrt(var_te)  > 0){r = 1}
          else{r=0}
        }
      }
      )
      unlist(power) |> na.omit() |> mean()
    }
    N <- c(25,50,100,150,200,500,1000,2000)#c(50,75,100,200,250,500)
    K<-6
    
    
    
    
    
    power_by_n <- unlist(sapply(N,function(n)calculate_power(n,TRUE)))

    power_Cluster <- data.frame(N=N,power_Cluster=power_by_n)
    
    #power_by_n <- unlist(sapply(N,function(n)calculate_power(n,TRUE))) #sfLapply

    #power_Cluster <- data.frame(N=N,power_Cluster=power_by_n)
    
    #power_by_n <- unlist(sapply(N, function(n)calculate_power(n,FALSE)))
    
    #power_HC <- data.frame(N=N,power_HC=power_by_n)
    #merge(power_Cluster,power_HC,by="N")
  })
  output$plot <- renderPlot({
    make_plot<-function(data){
      ggplot(data=data)+
        geom_line(mapping=aes(x=N,y=power_Cluster))+
        #geom_line(mapping=aes(x=N,y=power_HC),linetype = "dashed")+
        geom_hline(yintercept = 0.8)+
        scale_y_continuous(expand = c(0,0),,limits=c(0.,1),breaks = c(0,0.25,0.5,0.8))+
        scale_x_continuous(expand = c(0,0),limits=c(0,2010),breaks = c(0,100,250,500,1000,2000))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
        ylab(latex2exp::TeX("Power")) +
        xlab("Number of units (N)")
    }
    power <- powerVals()
    #plot(power$power_Cluster)
    make_plot(power)
  })
  
}

shinyApp(ui = ui, server = server)