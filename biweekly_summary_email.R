library(tidyverse)
library(rvest)
library(flextable)
library(officer)
library(RDCOMClient)  # Load the RDCOMClient library
library(knitr)
library(base64enc)
library(magick)

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
email[["To"]] <- "jaisrael@usbr.gov; cehlo@usbr.gov; avaisvil@usbr.gov; bmahardja@usbr.gov; lejohnson@usbr.gov; mmanzo@usbr.gov; dmmooney@usbr.gov; twashburn@usbr.gov"  # Change to the recipient's email address
email[["Subject"]] <- paste0("Salmonid Loss as of ", format(Sys.Date(), "%B %d, %Y"))

img_path <- file.path(tempdir(), "Loss_Table.png")
raster_img <- flextable::gen_grob(weekly_table)

# Save as PNG
png(img_path, width = 1000, height = 1200)
grid::grid.draw(raster_img)
dev.off()

# Attach to email
email[["Attachments"]]$Add(img_path)

# Save the plot as an image
img_path <- tempfile(fileext = ".png")
ggsave(img_path, plot = combined_graph, width = 10, height = 12)

# Read the image and convert to base64
img_base64 <- base64enc::dataURI(file = img_path, mime = "image/png")

# Construct the email body with embedded table and image
email_body <- paste0(
  "<p>Hi all,</p>",
  "<p>Please see the summary of the most recent loss estimates at salvage facilities.</p>",
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
  "<p>Season totals:</p>",
   "<ul>",
  "<li>Genetically confirmed: <strong>", final_summary$overall_summary$DNA_Run_W, "</strong></li>",
  "<li>Hatchery: <strong>", final_summary$overall_summary$CWT_Run_W, "</strong></li>",
  "<li>Unconfirmed LAD: <strong>", final_summary$overall_summary$unconfirmed_W, "</strong></li>",
  "</ul>",
  "<p>Previous week:</p>",
  "<ul>",
  "<li>Genetically confirmed: <strong>", final_summary$last_week_summary$DNA_Run_W, "</strong></li>",
  "<li>Hatchery: <strong>", final_summary$last_week_summary$CWT_Run_W, "</strong></li>",
  "<li>Unconfirmed LAD: <strong>", final_summary$last_week_summary$unconfirmed_W, "</strong></li>",
  "</ul>",
  "<p><img src='", img_base64, "' alt='Loss Graph' style='width: 100%; max-width: 800px;'/></p>",
  "<p>Best regards,<br>SaMT Team</p>" 
)

# Set the email body as HTML
email[["HTMLBody"]] <- email_body

# Save the email as a draft
email$Save()