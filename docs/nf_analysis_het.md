---
output: html_document
editor_options:
  chunk_output_type: console
css: highlight.css
---






# Analysis V (<i>H<sub>O</sub></i>)

## Summary

The population recombination rate is estimated within the [**nextflow**](https://www.nextflow.io/) script `analysis_het.nf` (located under `$BASE_DIR/nf/analysis_het/`), which runs on the _all BP_ data set.
Below is an overview of the steps involved in the analysis.
(The <span style="color:#4DAF4A">green dot</span> indicates the genotype input, <span style="color:#E41A1C">red arrows</span> depict output that is exported for further use.)

<div style="max-width:500px; margin:auto;">
<!--html_preserve--><div id="htmlwidget-df893c683a9272a004ab" style="width:672px;height:480px;" class="girafe html-widget"></div>
<script type="application/json" data-for="htmlwidget-df893c683a9272a004ab">{"x":{"html":"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e\" viewBox=\"0 0 432.00 432.00\">\n  <g>\n    <defs>\n      <clipPath id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_1\">\n        <rect x=\"0.00\" y=\"0.00\" width=\"432.00\" height=\"432.00\"/>\n      <\/clipPath>\n    <\/defs>\n    <rect x=\"0.00\" y=\"0.00\" width=\"432.00\" height=\"432.00\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_1\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_1)\" fill=\"#FFFFFF\" fill-opacity=\"1\" stroke-width=\"0.75\" stroke=\"#FFFFFF\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n    <defs>\n      <clipPath id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2\">\n        <rect x=\"0.00\" y=\"0.00\" width=\"432.00\" height=\"432.00\"/>\n      <\/clipPath>\n    <\/defs>\n    <g clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\">\n      <text x=\"111.63\" y=\"325.35\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_2\" font-size=\"225.00pt\" font-weight=\"bold\" fill=\"#E0E0E0\" fill-opacity=\"1\" font-family=\"DejaVu Sans\">7<\/text>\n    <\/g>\n    <polyline points=\"406.73,25.28 406.41,25.60 405.75,26.26 405.09,26.92 404.43,27.58 403.77,28.25 403.11,28.91 402.45,29.57 401.79,30.23 401.13,30.89 400.47,31.56 399.80,32.22 399.14,32.88 398.48,33.54 397.82,34.21 397.16,34.87 396.50,35.53 395.84,36.19 395.18,36.85 394.52,37.52 393.85,38.18 393.19,38.84 392.53,39.50 391.87,40.17 391.21,40.83 390.55,41.49 389.89,42.15 389.23,42.81 388.57,43.48 387.91,44.14 387.24,44.80 386.58,45.46 385.92,46.13 385.26,46.79 384.60,47.45 383.94,48.11 383.28,48.77 382.62,49.44 381.96,50.10 381.30,50.76 380.63,51.42 379.97,52.09 379.31,52.75 378.65,53.41 377.99,54.07 377.33,54.73 376.67,55.40 376.01,56.06 375.35,56.72 374.69,57.38 374.02,58.05 373.36,58.71 372.70,59.37 372.04,60.03 371.38,60.69 370.72,61.36 370.06,62.02 369.40,62.68 368.74,63.34 368.07,64.01 367.41,64.67 366.75,65.33 366.09,65.99 365.43,66.66 364.77,67.32 364.11,67.98 363.45,68.64 362.79,69.30 362.13,69.97 361.46,70.63 360.80,71.29 360.14,71.95 359.48,72.62 358.82,73.28 358.16,73.94 357.50,74.60 356.84,75.26 356.18,75.93 355.52,76.59 354.85,77.25 354.19,77.91 353.53,78.58 352.87,79.24 352.55,79.56\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_3\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"none\" stroke-width=\"1.06698\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polygon points=\"353.07,77.63 352.55,79.56 354.48,79.04\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_4\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"1.06698\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polyline points=\"341.29,90.84 340.97,91.16 340.31,91.82 339.65,92.48 338.99,93.14 338.33,93.80 337.67,94.47 337.01,95.13 336.35,95.79 335.68,96.45 335.02,97.11 334.36,97.78 333.70,98.44 333.04,99.10 332.38,99.76 331.72,100.43 331.06,101.09 330.40,101.75 329.73,102.41 329.07,103.07 328.41,103.74 327.75,104.40 327.09,105.06 326.43,105.72 325.77,106.38 325.11,107.05 324.45,107.71 323.79,108.37 323.12,109.03 322.46,109.69 321.80,110.36 321.14,111.02 320.48,111.68 319.82,112.34 319.16,113.00 318.50,113.67 317.84,114.33 317.17,114.99 316.51,115.65 315.85,116.31 315.19,116.98 314.53,117.64 313.87,118.30 313.21,118.96 312.55,119.62 311.89,120.29 311.23,120.95 310.56,121.61 309.90,122.27 309.24,122.93 308.58,123.60 307.92,124.26 307.26,124.92 306.60,125.58 305.94,126.25 305.28,126.91 304.62,127.57 303.95,128.23 303.29,128.89 302.63,129.56 301.97,130.22 301.31,130.88 300.65,131.54 299.99,132.20 299.33,132.87 298.67,133.53 298.00,134.19 297.34,134.85 296.68,135.51 296.02,136.18 295.36,136.84 294.70,137.50 294.04,138.16 293.38,138.82 292.72,139.49 292.06,140.15 291.39,140.81 290.73,141.47 290.07,142.13 289.41,142.80 288.75,143.46 288.09,144.12 287.43,144.78 287.11,145.10\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_5\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"none\" stroke-width=\"1.06698\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polygon points=\"287.62,143.18 287.11,145.10 289.03,144.58\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_6\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"1.06698\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polyline points=\"275.85,156.38 275.53,156.70 274.87,157.36 274.21,158.03 273.55,158.69 272.89,159.35 272.22,160.01 271.56,160.67 270.90,161.34 270.24,162.00 269.58,162.66 268.92,163.32 268.26,163.99 267.60,164.65 266.94,165.31 266.27,165.97 265.61,166.63 264.95,167.30 264.29,167.96 263.63,168.62 262.97,169.28 262.31,169.95 261.65,170.61 260.99,171.27 260.33,171.93 259.66,172.59 259.00,173.26 258.34,173.92 257.68,174.58 257.02,175.24 256.36,175.91 255.70,176.57 255.04,177.23 254.38,177.89 253.72,178.55 253.05,179.22 252.39,179.88 251.73,180.54 251.07,181.20 250.41,181.86 249.75,182.53 249.09,183.19 248.43,183.85 247.77,184.51 247.10,185.18 246.44,185.84 245.78,186.50 245.12,187.16 244.46,187.82 243.80,188.49 243.14,189.15 242.48,189.81 241.82,190.47 241.16,191.14 240.49,191.80 239.83,192.46 239.17,193.12 238.51,193.78 237.85,194.45 237.19,195.11 236.53,195.77 235.87,196.43 235.21,197.10 234.55,197.76 233.88,198.42 233.22,199.08 232.56,199.74 231.90,200.41 231.24,201.07 230.58,201.73 229.92,202.39 229.26,203.06 228.60,203.72 227.94,204.38 227.27,205.04 226.61,205.70 225.95,206.37 225.29,207.03 224.63,207.69 223.97,208.35 223.31,209.02 222.65,209.68 221.99,210.34 221.67,210.66\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_7\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"none\" stroke-width=\"1.06698\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polygon points=\"222.18,208.73 221.67,210.66 223.59,210.14\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_8\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"1.06698\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polyline points=\"210.40,221.93 210.09,222.24 209.42,222.90 208.76,223.56 208.10,224.22 207.44,224.88 206.78,225.54 206.12,226.20 205.46,226.86 204.79,227.52 204.13,228.18 203.47,228.84 202.81,229.50 202.15,230.16 201.49,230.82 200.83,231.48 200.17,232.14 199.50,232.80 198.84,233.46 198.18,234.13 197.52,234.79 196.86,235.45 196.20,236.11 195.54,236.77 194.88,237.43 194.21,238.09 193.55,238.75 192.89,239.41 192.23,240.07 191.57,240.73 190.91,241.39 190.25,242.05 189.59,242.71 188.92,243.37 188.26,244.03 187.60,244.69 186.94,245.35 186.28,246.01 185.62,246.67 184.96,247.33 184.30,247.99 183.63,248.65 182.97,249.31 182.31,249.97 181.65,250.63 180.99,251.29 180.33,251.95 179.67,252.61 179.01,253.27 178.34,253.93 177.68,254.59 177.02,255.25 176.36,255.91 175.70,256.57 175.04,257.23 174.38,257.89 173.72,258.55 173.05,259.21 172.39,259.87 171.73,260.53 171.07,261.19 170.41,261.85 169.75,262.51 169.09,263.17 168.42,263.83 167.76,264.49 167.10,265.16 166.44,265.82 165.78,266.48 165.12,267.14 164.46,267.80 163.80,268.46 163.13,269.12 162.47,269.78 161.81,270.44 161.15,271.10 160.49,271.76 159.83,272.42 159.17,273.08 158.51,273.74 157.84,274.40 157.18,275.06 156.52,275.72 156.21,276.03\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_9\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"none\" stroke-width=\"1.06698\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polygon points=\"156.73,274.10 156.21,276.03 158.14,275.52\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_10\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"1.06698\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polyline points=\"144.93,287.29 144.62,287.60 143.96,288.26 143.30,288.92 142.64,289.58 141.97,290.24 141.31,290.91 140.65,291.57 139.99,292.23 139.33,292.89 138.67,293.55 138.01,294.21 137.35,294.87 136.68,295.53 136.02,296.19 135.36,296.85 134.70,297.51 134.04,298.17 133.38,298.83 132.72,299.49 132.06,300.15 131.39,300.81 130.73,301.47 130.07,302.13 129.41,302.79 128.75,303.45 128.09,304.11 127.43,304.77 126.77,305.43 126.10,306.09 125.44,306.75 124.78,307.41 124.12,308.07 123.46,308.73 122.80,309.39 122.14,310.05 121.48,310.72 120.81,311.38 120.15,312.04 119.49,312.70 118.83,313.36 118.17,314.02 117.51,314.68 116.85,315.34 116.19,316.00 115.52,316.66 114.86,317.32 114.20,317.98 113.54,318.64 112.88,319.30 112.22,319.96 111.56,320.62 110.90,321.28 110.23,321.94 109.57,322.60 108.91,323.26 108.25,323.92 107.59,324.58 106.93,325.24 106.27,325.90 105.61,326.56 104.94,327.22 104.28,327.88 103.62,328.54 102.96,329.20 102.30,329.86 101.64,330.53 100.98,331.19 100.31,331.85 99.65,332.51 98.99,333.17 98.33,333.83 97.67,334.49 97.01,335.15 96.35,335.81 95.69,336.47 95.02,337.13 94.36,337.79 93.70,338.45 93.04,339.11 92.38,339.77 91.72,340.43 91.06,341.09 90.75,341.40\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_11\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"none\" stroke-width=\"1.06698\" stroke=\"#E41A1C\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polygon points=\"91.26,339.48 90.75,341.40 92.67,340.89\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_12\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"#E41A1C\" fill-opacity=\"1\" stroke-width=\"1.06698\" stroke=\"#E41A1C\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polyline points=\"79.46,352.66 79.15,352.97 78.49,353.63 77.83,354.29 77.17,354.95 76.51,355.61 75.85,356.27 75.19,356.93 74.53,357.59 73.86,358.25 73.20,358.91 72.54,359.57 71.88,360.23 71.22,360.89 70.56,361.55 69.90,362.21 69.23,362.87 68.57,363.53 67.91,364.19 67.25,364.85 66.59,365.51 65.93,366.17 65.27,366.83 64.61,367.49 63.94,368.15 63.28,368.81 62.62,369.47 61.96,370.13 61.30,370.79 60.64,371.45 59.98,372.11 59.31,372.77 58.65,373.43 57.99,374.09 57.33,374.75 56.67,375.41 56.01,376.07 55.35,376.73 54.69,377.39 54.02,378.05 53.36,378.71 52.70,379.37 52.04,380.03 51.38,380.69 50.72,381.35 50.06,382.01 49.40,382.67 48.73,383.33 48.07,383.99 47.41,384.65 46.75,385.31 46.09,385.97 45.43,386.63 44.77,387.29 44.10,387.95 43.44,388.61 42.78,389.27 42.12,389.93 41.46,390.59 40.80,391.25 40.14,391.91 39.48,392.57 38.81,393.23 38.15,393.89 37.49,394.55 36.83,395.21 36.17,395.87 35.51,396.53 34.85,397.19 34.19,397.85 33.52,398.51 32.86,399.17 32.20,399.83 31.54,400.49 30.88,401.15 30.22,401.81 29.56,402.47 28.89,403.13 28.23,403.78 27.57,404.44 26.91,405.10 26.25,405.76 25.59,406.42 25.28,406.73\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_13\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"none\" stroke-width=\"1.06698\" stroke=\"#E41A1C\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polygon points=\"25.80,404.81 25.28,406.73 27.20,406.22\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_14\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"#E41A1C\" fill-opacity=\"1\" stroke-width=\"1.06698\" stroke=\"#E41A1C\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <g clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\">\n      <text transform=\"translate(366.70,81.09) rotate(-405)\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_15\" font-size=\"8.28pt\" font-family=\"DejaVu Sans\">vcf_by_ind<\/text>\n    <\/g>\n    <g clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\">\n      <text transform=\"translate(307.14,140.74) rotate(-405)\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_16\" font-size=\"8.28pt\" font-family=\"DejaVu Sans\">inds_ch<\/text>\n    <\/g>\n    <g clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\">\n      <text transform=\"translate(108.98,338.87) rotate(-405)\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_17\" font-size=\"8.28pt\" font-family=\"DejaVu Sans\">inds_out<\/text>\n    <\/g>\n    <g clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\">\n      <text transform=\"translate(44.82,402.91) rotate(-405)\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_18\" font-size=\"8.28pt\" font-family=\"DejaVu Sans\">win_out<\/text>\n    <\/g>\n    <circle cx=\"412.36\" cy=\"19.64\" r=\"3.47pt\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_19\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"none\" stroke-width=\"0.708661\" stroke=\"#4DAF4A\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n    <circle cx=\"346.92\" cy=\"85.20\" r=\"3.47pt\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_20\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"none\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n    <circle cx=\"281.48\" cy=\"150.74\" r=\"3.47pt\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_21\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"none\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n    <circle cx=\"216.04\" cy=\"216.30\" r=\"3.47pt\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_22\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"none\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n    <circle cx=\"150.57\" cy=\"281.66\" r=\"3.47pt\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_23\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"none\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n    <circle cx=\"85.11\" cy=\"347.03\" r=\"3.47pt\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_24\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"none\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n    <circle cx=\"19.64\" cy=\"412.36\" r=\"3.47pt\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_25\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"none\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n    <circle cx=\"412.36\" cy=\"19.64\" r=\"1.87pt\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_26\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"#4DAF4A\" fill-opacity=\"1\" stroke-width=\"0.708661\" stroke=\"#4DAF4A\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\" title=\"Channel.fromFilePairs\"/>\n    <circle cx=\"346.92\" cy=\"85.20\" r=\"1.87pt\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_27\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\" title=\"split_inds\"/>\n    <circle cx=\"281.48\" cy=\"150.74\" r=\"1.87pt\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_28\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\" title=\"splitCsv\"/>\n    <circle cx=\"216.04\" cy=\"216.30\" r=\"1.87pt\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_29\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\" title=\"map\"/>\n    <circle cx=\"150.57\" cy=\"281.66\" r=\"1.87pt\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_30\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\" title=\"het_inds\"/>\n    <circle cx=\"85.11\" cy=\"347.03\" r=\"1.87pt\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_31\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\" title=\"win_inds\"/>\n    <circle cx=\"19.64\" cy=\"412.36\" r=\"1.87pt\" id=\"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_el_32\" clip-path=\"url(#svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\" title=\"\"/>\n  <\/g>\n<\/svg>\n","js":null,"uid":"svg_0c6ce85b-56e1-405e-9212-c7ea3f876d2e","ratio":1,"settings":{"tooltip":{"css":" .tooltip_SVGID_ { padding:5px;background:black;color:white;border-radius:2px 2px 2px 2px ; position:absolute;pointer-events:none;z-index:999;}\n","offx":10,"offy":0,"use_cursor_pos":true,"opacity":0.9,"usefill":false,"usestroke":false,"delay":{"over":200,"out":500}},"hover":{"css":" .hover_SVGID_ { fill:orange;stroke:gray; }\n"},"hoverkey":{"css":" .hover_key_SVGID_ { stroke:red; }\n"},"hovertheme":{"css":" .hover_theme_SVGID_ { fill:green; }\n"},"zoom":{"min":1,"max":1},"capture":{"css":" .selected_SVGID_ { fill:red;stroke:gray; }\n","type":"multiple","only_shiny":true,"selected":[]},"capturekey":{"css":" .selected_key_SVGID_ { stroke:gray; }\n","type":"single","only_shiny":true,"selected":[]},"capturetheme":{"css":" .selected_theme_SVGID_ { stroke:gray; }\n","type":"single","only_shiny":true,"selected":[]},"toolbar":{"position":"topright","saveaspng":true},"sizing":{"rescale":true,"width":1}}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->
</div>

## Details of `analysis_het.nf`

### Data preparation

The nextflow script starts by opening the genotype data.
Here we actually just do this to get a inventory of all genotyped samples - because of this we use the phased version (because it is the smaller file with the same samples) even though we are going to use the _all BP_ data set to actually estimate the heterozygosity.

<div class="kclass">

<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode">#<span class="hl opt">!/</span>usr<span class="hl opt">/</span>bin<span class="hl opt">/</span>env nextflow
<span class="hl slc">// nextflow run het.nf -c nextflow.config</span>
<span class="hl slc">// git 7.1</span>
<span class="hl slc">// open genotype data</span>
<span class="hl kwa">Channel</span>
	.fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/4_phased/phased_mac2.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> vcf_by_ind<span class="hl kwe">;</span> <span class="hl opt">}</span>
</code>
</pre>
</div>

The sample names are extracted from the genotype file and are fed into a channel.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 7.2</span>
<span class="hl slc">// collect all sample names</span>
<span class="hl kwa">process</span> split_inds <span class="hl opt">{</span>
	label <span class="hl str">&quot;L_loc_split_vcf&quot;</span>

	<span class="hl kwb">input</span><span class="hl kwe">:</span>
	<span class="hl kwd">set</span> vcfId<span class="hl opt">,</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">)</span> <span class="hl kwd">from</span> vcf_by_ind

	<span class="hl kwb">output</span><span class="hl kwe">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;inds.txt&quot;</span> <span class="hl opt">)</span> <span class="hl kwd">into</span> inds_ch

	<span class="hl kwb">script</span><span class="hl kwe">:</span>
	<span class="hl str">&quot;&quot;</span><span class="hl str">&quot;</span>
<span class="hl str">	vcfsamplenames ${vcf[0]} &gt; inds.txt</span>
<span class="hl str">	&quot;</span><span class="hl str">&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

The _all BP_ data set is split by sample, the allele counts for the sample are extracted and converted into a custom format.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 7.3</span>
<span class="hl slc">// for each sample, record allele frequencies</span>
<span class="hl kwa">process</span> het_inds <span class="hl opt">{</span>
	publishDir <span class="hl str">&quot;../../2_analysis/het/raw&quot;</span><span class="hl opt">,</span> mode<span class="hl kwe">:</span> <span class="hl str">&apos;copy&apos;</span>
	label <span class="hl str">&quot;L_190g15h_inds&quot;</span>

	<span class="hl kwb">input</span><span class="hl kwe">:</span>
	<span class="hl kwc">val</span><span class="hl opt">(</span> ind <span class="hl opt">)</span> <span class="hl kwd">from</span> inds_ch.splitCsv<span class="hl opt">()</span>.map<span class="hl opt">{</span><span class="hl kwc">it</span><span class="hl opt">[</span><span class="hl num">0</span><span class="hl opt">]}</span>

	<span class="hl kwb">output</span><span class="hl kwe">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*.het.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwd">into</span> inds_out

	<span class="hl kwb">script</span><span class="hl kwe">:</span>
	<span class="hl str">&quot;&quot;</span><span class="hl str">&quot;</span>
<span class="hl str">	vcftools --gzvcf \$BASE_DIR/1_genotyping/3_gatk_filtered/filterd_bi-allelic.allBP.vcf.gz \</span>
<span class="hl str">		--indv ${ind} \</span>
<span class="hl str">		--counts2 \</span>
<span class="hl str">		--stdout | \</span>
<span class="hl str">		sed  &apos;s/{COUNT}/REF</span><span class="hl esc">\\</span><span class="hl str">tALT/&apos; | \</span>
<span class="hl str">		cut -f 1,2,5,6 | \</span>
<span class="hl str">		awk -v OFS=&apos;</span><span class="hl esc">\\</span><span class="hl str">t&apos; -v ind=${ind} &apos;{if(NR==1){print \$1,\$2,&quot;HET&quot;,&quot;IND&quot;}else{x = \$3%2; ; print \$1,\$2,x/2,ind}}&apos; | \</span>
<span class="hl str">		gzip &gt; ${ind}.het.gz</span>
<span class="hl str">	&quot;</span><span class="hl str">&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

The heterozygosity is computed along non-operlapping 50 kb windows.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 7.3</span>
<span class="hl slc">// average heterozygosity over 50 kb windows</span>
<span class="hl kwa">process</span> win_inds <span class="hl opt">{</span>
	publishDir <span class="hl str">&quot;../../2_analysis/het/50kb&quot;</span><span class="hl opt">,</span> mode<span class="hl kwe">:</span> <span class="hl str">&apos;copy&apos;</span>
	label <span class="hl str">&quot;L_20g2h_inds&quot;</span>
	module <span class="hl str">&quot;R3.5.2&quot;</span>

	<span class="hl kwb">input</span><span class="hl kwe">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> het <span class="hl opt">)</span> <span class="hl kwd">from</span> inds_out

	<span class="hl kwb">output</span><span class="hl kwe">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;*_win_het.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwd">into</span> win_out

	<span class="hl kwb">script</span><span class="hl kwe">:</span>
	<span class="hl str">&quot;&quot;</span><span class="hl str">&quot;</span>
<span class="hl str">	Rscript --vanilla \$BASE_DIR/R/het_by_ind.R ${het} 50000</span>
<span class="hl str">	&quot;</span><span class="hl str">&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>
</div>

Finally, we are done with heterozygosity.

---
