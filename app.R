library(shiny)
library(plotly)


ui <- navbarPage(title = "ShinyGPAS",
                 tabPanel(title = "S1",
                          sidebarPanel(
                            sliderInput("h21",
                                        "Heritability (h2):",
                                        min = 0,
                                        max = 1,
                                        value = 0.5),
                            
                            sliderInput("N1",
                                        "Number of individuals (N):",
                                        min = 0,
                                        max = 30000,
                                        value = 1000),
                            
                            sliderInput("Me1",
                                        "Number of independent chromosome segments (Me):",
                                        min = 0,
                                        max = 10000,
                                        value = 1000),
                            
                            radioButtons('xaxis', 'X-axis', c("h2"="h21", "N"="N1", "Me"="Me1"), inline = TRUE)
                          ),
                          
                          tags$a(href="https://doi.org/10.1371/journal.pone.0003395", "[10.1371/journal.pone.0003395]"), tags$a(href="https://doi.org/10.1534/genetics.110.116855","[10.1534/genetics.110.116855]"),
                          
                          mainPanel(
                            plotlyOutput("distPlot")
                          )
                 ),
                 
                 tabPanel(title = "S2",
                          sidebarPanel(
                            sliderInput("h22",
                                        "Heritability (h2):",
                                        min = 0,
                                        max = 1,
                                        value = 0.5),
                            
                            sliderInput("N2",
                                        "Number of individuals (N):",
                                        min = 0,
                                        max = 30000,
                                        value = 1000),
                            
                            sliderInput("Me2",
                                        "Number of independent chromosome segments (Me):",
                                        min = 0,
                                        max = 10000,
                                        value = 1000),
                            
                            sliderInput("Ne2",
                                        "Number of effective population size (Ne):",
                                        min = 0,
                                        max = 10000,
                                        value = 500),
                            
                            radioButtons('xaxis2', 'X-axis', c("h2"="h22", "N"="N2", "Me" = "Me2", "Ne" = "Ne2"), inline = TRUE)
                          ),
                          
                          tags$a(href="https://doi.org/10.1007/s10709-008-9308-0", "[10.1007/s10709-008-9308-0]"), 
                          
                          mainPanel(
                            plotlyOutput("distPlot2")
                          )
                 ),
                 
                 tabPanel(title = "S3",
                          sidebarPanel(
                            sliderInput("h23",
                                        "Heritability (h2):",
                                        min = 0,
                                        max = 1,
                                        value = 0.5),
                            
                            sliderInput("N3",
                                        "Number of individuals (N):",
                                        min = 0,
                                        max = 30000,
                                        value = 1000),
                            
                            sliderInput("Me3",
                                        "Number of independent chromosome segments (Me):",
                                        min = 0,
                                        max = 10000,
                                        value = 1000),
                            
                            sliderInput("b3",
                                        "Proportion of genetic variance explained by the markers (b):",
                                        min = 0,
                                        max = 1,
                                        value = 0.5),
                            
                            radioButtons('xaxis3', 'X-axis', c("h2"="h23", "N"="N3", "Me" = "Me3", "b" = "b3"), inline = TRUE)
                          ),
                          
                          tags$a(href="https://doi.org/10.1111/j.1439-0388.2011.00964.x", "[10.1111/j.1439-0388.2011.00964.x]"), 
                          mainPanel(
                            plotlyOutput("distPlot3")
                          )
                 ),
                 

                 
                 tabPanel(title = "S4",
                          sidebarPanel(
                            sliderInput("h25",
                                        "Heritability (h2):",
                                        min = 0,
                                        max = 1,
                                        value = 0.5),
                            
                            sliderInput("b5",
                                        "Regression of genomic relationships realized at markers on genetic relationships realized at causal loci (b):",
                                        min = 0,
                                        max = 1,
                                        value = 0.5),
                            
                            radioButtons('xaxis5', 'X-axis', c("h2"="h25", "b" = "b5"), inline = TRUE)
                          ),
                          
                          tags$a(href="https://doi.org/10.1371/journal.pgen.1003608", "[10.1371/journal.pgen.1003608]"),
                          
                          mainPanel(
                            plotlyOutput("distPlot5")
                          )
                 ),
                 
                 
                 tabPanel(title = "S5",
                          sidebarPanel(
                            sliderInput("h24",
                                        "Genomic Heritability (h2):",
                                        min = 0,
                                        max = 1,
                                        value = 0.5),
                            
                            sliderInput("N4",
                                        "Number of individuals (N):",
                                        min = 0,
                                        max = 30000,
                                        value = 1000),
                            
                            sliderInput("Me4",
                                        "Number of independent chromosome segments (Me):",
                                        min = 0,
                                        max = 10000,
                                        value = 1000),
                            
                            radioButtons('xaxis4', 'X-axis', c("h2"="h24", "N"="N4", "Me"="Me4"), inline = TRUE)
                          ),
                          
                          tags$a(href="https://doi.org/10.1371/journal.pone.0161054", "[10.1371/journal.pone.0161054]"),
                          
                          mainPanel(
                            plotlyOutput("distPlot4")
                          )
                 ),
                 
                 
                 
                 
                 
                 tabPanel(title = "S6",
                          sidebarPanel(
                            
                            sliderInput("rG_AB6",
                                        "Genetic correlation between populations A and B (rG_AB):",
                                        min = 0,
                                        max = 1,
                                        value = 0.5),
                            
                            
                            sliderInput("rG_AC6",
                                        "Genetic correlation between populations A and C (rG_AC):",
                                        min = 0,
                                        max = 1,
                                        value = 0.5),
                            
                            sliderInput("rG_BC6",
                                        "Genetic correlation between populations B and C (rG_BC):",
                                        min = 0,
                                        max = 1,
                                        value = 0.5),
                            
                            sliderInput("b_AC6",
                                        "Square root of the proportion of the genetic variance in predicted population C explained by the markers in population A (b_AC):",
                                        min = 0,
                                        max = 1,
                                        value = 0.5),
                            
                            sliderInput("b_BC6",
                                        "Square root of the proportion of the genetic variance in predicted population C explained by the markers in population B (b_BC):",
                                        min = 0,
                                        max = 1,
                                        value = 0.5),
                            
                           
                            
                            
                           
                            
                            sliderInput("h26a",
                                        "Heritability in population A (h2):",
                                        min = 0,
                                        max = 1,
                                        value = 0.5),
                            
                            sliderInput("h26b",
                                        "Heritability in population B (h2):",
                                        min = 0,
                                        max = 1,
                                        value = 0.5),
                            
                            sliderInput("N6a",
                                        "Number of individuals in population A (N):",
                                        min = 0,
                                        max = 30000,
                                        value = 1000),
                            
                            sliderInput("N6b",
                                        "Number of individuals in population B(N):",
                                        min = 0,
                                        max = 30000,
                                        value = 1000),
                            
                            
                            sliderInput("Me6",
                                        "Number of independent chromosome segments (Me):",
                                        min = 0,
                                        max = 10000,
                                        value = 1000),
                            
                            
                            radioButtons('xaxis6', 'X-axis', c( "rG" = "rG_AC6", "b" = "b_AC6",  "h2"="h26a", "N" = "N6a"), inline = TRUE)
                          ),
                          
                          tags$a(href="https://doi.org/10.1534/genetics.115.183269", "[10.1534/genetics.115.1832698]"),
                          
                          mainPanel(
                            plotlyOutput("distPlot6")
                          )
                 )
                 
                 
)


server <- function(input, output) {
  
    # S1 
    output$distPlot <- renderPlotly({
    if (input$xaxis == "h21") {
    h2 <- seq(0, 1, by=0.001)
    r <- sqrt(  (input$N1 * h2)/(input$N1 * h2 + input$Me1)  ) 
    dat <- data.frame(r=r, h2=h2)
    return(plot_ly(dat, x=~h2, y=~r, type = 'scatter', mode = 'markers', text = ~paste('h2:', h2, "<br>r: ", r)))
    } else if (input$xaxis == "N1") {
    N <- 0:input$N1
    r <- sqrt(  (N * input$h21)/(N * input$h21 + input$Me1)  ) 
    dat <- data.frame(r=r, N=N)
    return(plot_ly(dat, x=~N, y=~r, type = 'scatter', mode = 'markers',text = ~paste('N:', N, "<br>r: ", r))) 
    } else (input$xaxis == "Me1") 
    Me <- 0:input$Me1
    r <- sqrt(  (input$N1 * input$h21)/(input$N1  * input$h21 + Me)  ) 
    dat <- data.frame(r=r, Me=Me)
    return(plot_ly(dat, x=~Me , y=~r, type = 'scatter', mode = 'markers',text = ~paste('Me:', Me, "<br>r: ", r)))
    
  })

    
    # S2
    output$distPlot2 <- renderPlotly({
    if (input$xaxis2 == "h22") {
    h2 <- seq(0, 1, by=0.001)
    lambda <- input$Me2 / (h2 * log(2 * input$Ne2))
    alpha <- 1 + 2 *( input$Me2 / (input$N2 * h2 * log(2 * input$Ne2)) )
    r <- sqrt(  1 - lambda/(2 * input$N2 * sqrt(alpha)) * log( (1 + alpha + 2 * sqrt(alpha))/(1 + alpha - 2 * sqrt(alpha)) )  ) 
    dat <- data.frame(r=r, h2=h2)
    return(plot_ly(dat, x=~h2, y=~r, type = 'scatter', mode = 'markers', text = ~paste('h2:', h2, "<br>r: ", r)))
    } else if (input$xaxis2 == "N2") {
    N <- 0:input$N2
    lambda <- input$Me2 / (input$h22 * log(2 * input$Ne2))
    alpha <- 1 + 2 *( input$Me2 / (N * input$h22 * log(2 * input$Ne2)) )
    r <- sqrt(  1 - lambda/(2 * N * sqrt(alpha)) * log( (1 + alpha + 2 * sqrt(alpha))/(1 + alpha - 2 * sqrt(alpha)) )  ) 
    dat <- data.frame(r=r, N=N)
    return(plot_ly(dat, x=~N, y=~r, type = 'scatter', mode = 'markers',text = ~paste('N:', N, "<br>r: ", r))) 
    } else if (input$xaxis2 == "Me2"){ 
    Me <- 0:input$Me2
    lambda <- Me / (input$h22 * log(2 * input$Ne2))
    alpha <- 1 + 2 *(Me / (input$N2 * input$h22 * log(2 * input$Ne2)) )
    r <- sqrt(  1 - lambda/(2 * input$N2 * sqrt(alpha)) * log( (1 + alpha + 2 * sqrt(alpha))/(1 + alpha - 2 * sqrt(alpha)) )  ) 
    dat <- data.frame(r=r, Me=Me)
    return(plot_ly(dat, x=~Me , y=~r, type = 'scatter', mode = 'markers',text = ~paste('Me:', Me, "<br>r: ", r)))
    } else  (input$xaxis2 == "Ne2") 
    Ne <- 0:input$Ne2
    lambda <- input$Me2 / (input$h22 * log(2 * Ne))
    alpha <- 1 + 2 *(input$Me2 / (input$N2 * input$h22 * log(2 * Ne)) )
    r <- sqrt(  1 - lambda/(2 * input$N2 * sqrt(alpha)) * log( (1 + alpha + 2 * sqrt(alpha))/(1 + alpha - 2 * sqrt(alpha)) )  ) 
    dat <- data.frame(r=r, Ne=Ne)
    return(plot_ly(dat, x=~Ne , y=~r, type = 'scatter', mode = 'markers',text = ~paste('Ne:', Ne, "<br>r: ", r)))
      
    })
    
    # S3
    output$distPlot3 <- renderPlotly({
    if (input$xaxis3 == "h23") {
    h2 <- seq(0, 1, by=0.001)
    tmp3 <- (input$N3 * input$b3 * h2) /input$Me3
    r <- sqrt( input$b3 *  ((tmp3) / (1 + tmp3)) ) 
    dat <- data.frame(r=r, h2=h2)
    return(plot_ly(dat, x=~h2, y=~r, type = 'scatter', mode = 'markers', text = ~paste('h2:', h2, "<br>r: ", r)))
    } else if (input$xaxis3 == "N3") {
    N <- 0:input$N3
    tmp3 <- (N * input$b3 * input$h23) /input$Me3
    r <- sqrt( input$b3 *  ((tmp3) / (1 + tmp3)) ) 
    dat <- data.frame(r=r, N=N)
    return(plot_ly(dat, x=~N, y=~r, type = 'scatter', mode = 'markers',text = ~paste('N:', N, "<br>r: ", r))) 
    } else if (input$xaxis3 == "Me3"){ 
    Me <- 0:input$Me3
    tmp3 <- (input$N3 * input$b3 * input$h23) /Me
    r <- sqrt( input$b3 *  ((tmp3) / (1 + tmp3)) ) 
    dat <- data.frame(r=r, Me=Me)
    return(plot_ly(dat, x=~Me , y=~r, type = 'scatter', mode = 'markers',text = ~paste('Me:', Me, "<br>r: ", r)))
    } else  (input$xaxis3 == "b3") 
    b <- seq(0, 1, by=0.001)
    tmp3 <- (input$N3 * b * input$h23) / input$Me3
    r <- sqrt( b *  ((tmp3) / (1 + tmp3)) ) 
    dat <- data.frame(r=r, b=b)
    return(plot_ly(dat, x=~b , y=~r, type = 'scatter', mode = 'markers',text = ~paste('b:', b, "<br>r: ", r)))
        
  })
    
  
    
    
    # S5
    output$distPlot5 <- renderPlotly({
    if (input$xaxis5 == "h25") {
    h2 <- seq(0, 1, by=0.001)
    R2 <-  (1 - (1 - input$b5)^2 ) * h2    
    dat <- data.frame(R2=R2, h2=h2)
    return(plot_ly(dat, x=~h2, y=~R2, type = 'scatter', mode = 'markers', text = ~paste('h2:', h2, "<br>R2: ", R2)))
    }  else  (input$xaxis5 == "b5") 
    b <- seq(0, 1, by=0.001)
    R2 <-  (1 - (1 - b)^2 ) * input$h25   
    dat <- data.frame(R2=R2, b=b)
    return(plot_ly(dat, x=~b , y=~R2, type = 'scatter', mode = 'markers',text = ~paste('b:', b, "<br>R2: ", R2)))
        
  })
    
    
    # S4 
    output$distPlot4 <- renderPlotly({
      if (input$xaxis4 == "h24") {
        h2 <- seq(0, 1, by=0.001)
        R2 <-   h2 * ((input$N4 * h2)/(input$N4 * h2 + input$Me4))   
        dat <- data.frame(R2=R2, h2=h2)
        return(plot_ly(dat, x=~h2, y=~R2, type = 'scatter', mode = 'markers', text = ~paste('h2:', h2, "<br>R2: ", R2)))
      } else if (input$xaxis4 == "N4") {
        N <- 0:input$N4
        R2 <-   input$h24 * ((N * input$h24)/(N * input$h24 + input$Me4))
        dat <- data.frame(R2=R2, N=N)
        return(plot_ly(dat, x=~N, y=~R2, type = 'scatter', mode = 'markers',text = ~paste('N:', N, "<br>R2: ", R2))) 
      } else (input$xaxis4 == "Me4") 
      Me <- 0:input$Me4
      R2 <-  input$h24 * ((input$N4 * input$h24)/(input$N4  * input$h24 + Me) )
      dat <- data.frame(R2=R2, Me=Me)
      return(plot_ly(dat, x=~Me , y=~R2, type = 'scatter', mode = 'markers',text = ~paste('Me:', Me, "<br>R2: ", R2)))
      
    })
    
    # S6
    output$distPlot6 <- renderPlotly({
       if (input$xaxis6 == "b_AC6") {
         b_AC <- seq(0, 1, by=0.1)
         b_BC <- seq(0, 1, by=0.1)
        
        fun2 <- function(b_AC, b_BC){
          
          term1 <- matrix(c(b_AC * input$rG_AC6 * sqrt(input$h26a/input$Me6), b_BC * input$rG_BC6 * sqrt(input$h26b/input$Me6)), nrow=1, ncol=2)
          
          term2 <- matrix(c(input$h26a/input$Me6 + 1/input$N6a, input$rG_AB6 * sqrt(input$h26a*input$h26b)/input$Me6,  input$rG_AB6 * sqrt(input$h26a*input$h26b)/input$Me6, input$h26b/input$Me6 + 1/input$N6b ), nrow=2, ncol=2) 
          term2 <- term2 + diag(ncol(term2))*0.001 
          term2 <- solve(term2)
          
          term3 <- matrix(c(b_AC * input$rG_AC6 * sqrt(input$h26a/input$Me6), b_BC * input$rG_BC6 * sqrt(input$h26b/input$Me6)), nrow=2, ncol=1)
          
          ri <- sqrt( term1 %*% term2 %*% term3  ) 
          return(ri)
        }
        
        r <- array()
        b_AC.vec <- array()
        b_BC.vec <- array()
        index <- 1
        for (i in 1:length(b_AC)){
          for (j in 1:length(b_BC)){
            b_AC.vec[index] <- b_AC[i]
            b_BC.vec[index] <- b_BC[j]
            r[index] <- fun2(b_AC=b_AC[i], b_BC=b_BC[j])
            index <- index + 1
          }
        }
        
        dat <- data.frame(r=r, b_AC=b_AC.vec, b_BC=b_BC.vec)
        
        return(plot_ly(dat, x=~b_AC, y=~b_BC, z=~r, type="scatter3d", text = ~paste('b_AC:', b_AC, "<br>b_AC: ", b_BC, "<br>r: ", r))) 
       }  else  if (input$xaxis6 == "rG_AC6"){ 
      rG_AC <- seq(0, 1, by=0.1)
      rG_BC <- seq(0, 1, by=0.1)
      
      fun2 <- function(rG_AC, rG_BC){
        
        term1 <- matrix(c(input$b_AC6 * rG_AC * sqrt(input$h26a/input$Me6), input$b_BC6 * rG_BC * sqrt(input$h26b/input$Me6)), nrow=1, ncol=2)
        
        term2 <- matrix(c(input$h26a/input$Me6 + 1/input$N6a, input$rG_AB6 * sqrt(input$h26a*input$h26b)/input$Me6,  input$rG_AB6 * sqrt(input$h26a*input$h26b)/input$Me6, input$h26b/input$Me6 + 1/input$N6b ), nrow=2, ncol=2) 
        term2 <- term2 + diag(ncol(term2))*0.001 
        term2 <- solve(term2)
        
        term3 <- matrix(c(input$b_AC6 * rG_AC * sqrt(input$h26a/input$Me6), input$b_BC6 * rG_BC * sqrt(input$h26b/input$Me6)), nrow=2, ncol=1)
        
        ri <- sqrt( term1 %*% term2 %*% term3  ) 
        return(ri)
      }
      
      r <- array()
      rG_AC.vec <- array()
      rG_BC.vec <- array()
      index <- 1
      for (i in 1:length(rG_AC)){
        for (j in 1:length(rG_BC)){
          rG_AC.vec[index] <- rG_AC[i]
          rG_BC.vec[index] <- rG_BC[j]
          r[index] <- fun2(rG_AC=rG_AC[i], rG_BC=rG_BC[j])
          index <- index + 1
        }
      }
      
      dat <- data.frame(r=r, rG_AC=rG_AC.vec, rG_BC=rG_BC.vec)
      
      return(plot_ly(dat, x=~rG_AC, y=~rG_BC, z=~r, type="scatter3d", text = ~paste('rG_AC:', rG_AC, "<br>rG_BC: ", rG_BC, "<br>r: ", r))) 
      
    } else if (input$xaxis6 == "h26a") {
        h2a <- seq(0, 1, by=0.1)
        h2b <- seq(0, 1, by=0.1)
        
        fun1 <- function(h2a, h2b){
          
          term1 <- matrix(c(input$b_AC6 * input$rG_AC6 * sqrt(h2a/input$Me6), input$b_BC6 * input$rG_BC6 * sqrt(h2b/input$Me6)), nrow=1, ncol=2)
          
          term2 <- matrix(c(h2a/input$Me6 + 1/input$N6a, input$rG_AB6 * sqrt(h2a*h2b)/input$Me6,  input$rG_AB6 * sqrt(h2a*h2b)/input$Me6, h2b/input$Me6 + 1/input$N6b ), nrow=2, ncol=2) 
          term2 <- term2 + diag(ncol(term2))*0.001 
          term2 <- solve(term2)
          
          term3 <- matrix(c(input$b_AC6 * input$rG_AC6 * sqrt(h2a/input$Me6), input$b_BC6 * input$rG_BC6 * sqrt(h2b/input$Me6)), nrow=2, ncol=1)
          
          ri <- sqrt( term1 %*% term2 %*% term3  ) 
          return(ri)
        }
        
        r <- array()
        h2a.vec <- array()
        h2b.vec <- array()
        index <- 1
        for (i in 1:length(h2a)){
          for (j in 1:length(h2b)){
          h2a.vec[index] <- h2a[i]
          h2b.vec[index] <- h2b[j]
          r[index] <- fun1(h2a=h2a[i], h2b=h2b[j])
          index <- index + 1
          }
        }
        
        dat <- data.frame(r=r, h2a=h2a.vec, h2b=h2b.vec)
        
        return(plot_ly(dat, x=~h2a, y=~h2b, z=~r, type="scatter3d", text = ~paste('h2a:', h2a, "<br>h2b: ", h2b, "<br>r: ", r))) 

      }  else (input$xaxis6 == "N6a") 
        Na <- seq(1, 5000, by = 500)
        Nb <- seq(1, 5000, by = 500)
        
        fun2 <- function(Na, Nb){
          
          term1 <- matrix(c(input$b_AC6 * input$rG_AC6 * sqrt(input$h26a/input$Me6), input$b_BC6 * input$rG_BC6 * sqrt(input$h26b/input$Me6)), nrow=1, ncol=2)
          
          term2 <- matrix(c(input$h26a/input$Me6 + 1/Na, input$rG_AB6 * sqrt(input$h26a*input$h26b)/input$Me6,  input$rG_AB6 * sqrt(input$h26a*input$h26b)/input$Me6, input$h26b/input$Me6 + 1/Nb ), nrow=2, ncol=2) 
          term2 <- term2 + diag(ncol(term2))*0.001 
          term2 <- solve(term2)
          
          term3 <- matrix(c(input$b_AC6 * input$rG_AC6 * sqrt(input$h26a/input$Me6), input$b_BC6 * input$rG_BC6 * sqrt(input$h26b/input$Me6)), nrow=2, ncol=1)
          
          ri <- sqrt( term1 %*% term2 %*% term3  ) 
          return(ri)
        }
        
        r <- array()
        Na.vec <- array()
        Nb.vec <- array()
        index <- 1
        for (i in 1:length(Na)){
          for (j in 1:length(Nb)){
            Na.vec[index] <- Na[i]
            Nb.vec[index] <- Nb[j]
            r[index] <- fun2(Na=Na[i], Nb=Nb[j])
            index <- index + 1
          }
        }
        
        dat <- data.frame(r=r, Na=Na.vec, Nb=Nb.vec)
        
        return(plot_ly(dat, x=~Na, y=~Nb, z=~r, type="scatter3d", text = ~paste('Na:', Na, "<br>Nb: ", Nb, "<br>r: ", r))) 
      
        
        
        
    })
    
    
}



shinyApp(ui = ui, server = server)


