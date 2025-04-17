library(tidyverse)
library(rvest)
library(flextable)
library(officer)
library(RDCOMClient)  # Load the RDCOMClient library
library(knitr)
library(base64enc)
library(magick)
library(kableExtra)

source('LossSummary.R')

# Calculate loss estimates
WRloss <- cumulative_loss %>% filter(species == 'Winter-run') %>% summarize(maxloss = max(cumul_loss)) %>% pull()
partial <- temp %>% filter(confirmed == 'PARTIAL') %>% summarize(sum(loss)) %>% pull()
SHloss <- cumulative_loss %>% filter(species == 'Steelhead') %>% summarize(maxloss = max(cumul_loss)) %>% pull()
WRtriggers <- sum(wr_weekly$triggered == "Yes")
SHtriggers <- sum(SH_weekly$triggered == "Yes")
WR_hatch_loss <- max(wr_hatch$cumul_loss)

# Create a new Outlook application object
outlook_app <- COMCreate("Outlook.Application")

# Create a new mail item
email <- outlook_app$CreateItem(0)  # 0 corresponds to olMailItem

# Set email properties
email[["To"]] <- "jaisrael@usbr.gov; cehlo@usbr.gov; avaisvil@usbr.gov; bmahardja@usbr.gov; lejohnson@usbr.gov; kristin.begun@noaa.gov;
mmanzo@usbr.gov; dmmooney@usbr.gov; twashburn@usbr.gov; Farida.Islam@water.ca.gov; Jeffrey.Onsted@water.ca.gov; ashamilton@usbr.gov"  # Change to the recipient's email address
email[["Subject"]] <- paste0("Salmonid Loss as of ", format(Sys.Date(), "%B %d, %Y"))

# Save the plot as an image
img_path <- tempfile(fileext = ".png")
ggsave(img_path, plot = combined_graph, width = 10, height = 12)

# Read the image and convert to base64
img_base64 <- base64enc::dataURI(file = img_path, mime = "image/png")

n_cols <- ncol(SH_weekly)
column_width <- 100  # adjust to fit your content

sh_weekly_html <- knitr::kable(
  SH_weekly, format = "html",
  align = rep("l", n_cols),
  table.attr = "style='border-collapse: collapse; margin: 0; padding: 0; line-height: 1.2; width: 700px; font-family: Arial; font-size: 10pt;'"
) %>%
  kable_styling(
    full_width = FALSE,
    position = "left",
    font_size = 12,
    stripe_color = "#f9f9f9"
  ) %>%
  row_spec(0, bold = TRUE,
           extra_css = "border-top: 1px solid black; border-bottom: 1px solid black;") %>%
  column_spec(
    1:n_cols,
    width = paste0(column_width, "px"),
    extra_css = "border-left: none; border-right: none;"
  ) %>%
  # Add bottom border to last row
  row_spec(nrow(SH_weekly), extra_css = "border-bottom: 1px solid black;")

# Construct the email body with embedded SH_weekly data
email_body <- paste0(
  "<p>Hi all,</p>",
  "<p>Please see the summary of the most recent loss estimates at salvage facilities and please note that data is preliminary and subject to change.</p>",
  "<ul>",
  "<li>Total annual loss of natural Winter-run Chinook Salmon as of ", 
  format(Sys.Date() - 1, "%B %d, %Y"), 
  " is <strong>", WRloss - partial, "</strong>, and loss of unconfirmed Winter-run Chinook Salmon is <strong>", partial, "</strong>.</li>",
  "<li>Total annual loss of hatchery Winter-run Chinook Salmon as of ", 
  format(Sys.Date() - 1, "%B %d, %Y"), 
  " is <strong>", WR_hatch_loss, "</strong>.</li>",
  "<li>Total annual loss of CCV steelhead as of ", 
  format(Sys.Date() - 1, "%B %d, %Y"), 
  " is <strong>", SHloss, "</strong>.</li>",
  "<li>There have been <strong>", WRtriggers, "</strong> triggers of the Winter-run and <strong>", SHtriggers, "</strong> triggers of the steelhead distributed weekly loss thresholds in the past week.</li>",
  "<ul>",
  "<p>Season total Winter-run Chinook salmon:</p>",
  "<ul>",
  "<li>Genetically confirmed: <strong>", final_summary$overall_summary$DNA_Run_W, "</strong></li>",
  "<li>Hatchery: <strong>", final_summary$overall_summary$CWT_Run_W, "</strong></li>",
  "<li>Unconfirmed LAD: <strong>", final_summary$overall_summary$unconfirmed_W, "</strong></li>",
  "</ul>",
  "<p>Previous week Winter-run Chinook salmon:</p>",
  "<ul>",
  "<li>Genetically confirmed: <strong>", final_summary$last_week_summary$DNA_Run_W, "</strong></li>",
  "<li>Hatchery: <strong>", final_summary$last_week_summary$CWT_Run_W, "</strong></li>",
  "<li>Unconfirmed LAD: <strong>", final_summary$last_week_summary$unconfirmed_W, "</strong></li>",
  "</ul>",
  "<p>Weekly Steelhead Data:</p>",
  sh_weekly_html,
  "<p><img src='", img_base64, "' alt='Loss Graph' style='width: 100%; max-width: 800px;'/></p>",
  "<p>Best regards,<br>SaMT Team</p>" 
)

# Set the email body as HTML
email[["HTMLBody"]] <- email_body

# Save the email as a draft
email$Save()