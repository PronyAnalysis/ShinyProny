library(shiny)
library(ggplot2)

shinyUI(pageWithSidebar(
  
  headerPanel("Prony analysis"),
  
  sidebarPanel(    
    sliderInput('mu1', 'mu1', min = 0.01, max = 5, value = 1, step = 0.01),
    sliderInput('mu2', 'mu2', min = 0.01, max = 5, value = 1.5, step = 0.01),
    sliderInput('a1', 'a1', min = 0.01, max = 5, value = 1, step = 0.01),
    sliderInput('a2', 'a2', min = 0.01, max = 5, value = 1, step = 0.01),
    sliderInput('delta', 'delta', min = 0.01, max = 0.5, value = 0.05, step = 0.01),
    sliderInput('noise', 'noise', min = 0.00, max = 1, value = 0.1, step = 0.01)
  ),
  
  mainPanel(
    plotOutput('plot', height = 700),
    dataTableOutput("mytable"),
    helpText(a("Click here to view the detailed system description", href = "https://www.authorea.com/users/44760/articles/51974/_show_article"))
  )
)) 