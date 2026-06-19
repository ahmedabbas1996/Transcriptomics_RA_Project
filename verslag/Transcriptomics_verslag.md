# Transcriptomische analyse van reumatoïde artritis met RNA-sequencing

## Inhoudsopgave
- [Introductie](#introductie)
- [Methoden](#methoden)
- [Resultaten](#resultaten)
- [Conclusie](#conclusie)
- [Gebruik van GitHub](#gebruik-van-github)
- [Referenties](#referenties)

---

## Introductie

Reumatoïde artritis (RA) is een chronische auto-immuunziekte die wordt gekenmerkt door ontsteking van de gewrichten en geleidelijke beschadiging van kraakbeen en bot. De ziekte ontstaat door een ontregelde immuunrespons waarbij verschillende immuuncellen, cytokinen en signaalroutes betrokken zijn. Ondanks beschikbare behandelingen is nog niet volledig bekend welke moleculaire processen betrokken zijn bij het ontstaan en de progressie van RA (Smolen et al., 2018).

Met behulp van RNA-sequencing (RNA-seq) kan de expressie van duizenden genen tegelijkertijd worden onderzocht. Hierdoor kunnen verschillen tussen patiënten en gezonde controles worden geïdentificeerd. Vervolgens kunnen functionele analyses, zoals Gene Ontology (GO) and KEGG pathwayanalyse, worden gebruikt om de biologische betekenis van deze genexpressieveranderingen te interpreteren (Conesa et al., 2016).

Het doel van dit onderzoek was het identificeren van differentieel geëxprimeerde genen tussen RA-patiënten en gezonde controles. Daarnaast werd onderzocht welke biologische processen en signaalroutes verrijkt zijn in RA. De onderzoeksvraag luidde:

**Welke genexpressieveranderingen en biologische pathways zijn geassocieerd met reumatoïde artritis?**

---

## Methoden

Voor deze analyse werd gebruikgemaakt van een openbare RNA-seq dataset bestaande uit vier RA-monsters en vier gezonde controles ($n = 8$). De ruwe countmatrix werd ingelezen in R en gecombineerd met de bijbehorende samplemetadata.

De differentiële genexpressieanalyse werd uitgevoerd met het pakket `DESeq2`. Genen werden als significant beschouwd bij een aangepaste p-waarde ($padj < 0.05$). Voor de volcano plot werden genen met een absolute log2 fold change ($|\log_2\text{FC}| > 1$) extra benadrukt.

Om de globale verschillen tussen de monsters te onderzoeken werd een Principal Component Analysis (PCA) uitgevoerd. De resultaten van de differentiële expressieanalyse werden gevisualiseerd met een volcano plot.

Voor de biologische interpretatie van de resultaten werd een Gene Ontology (GO) Biological Process analyse uitgevoerd. Daarnaast werd een KEGG pathwayanalyse uitgevoerd met behulp van het pakket `clusterProfiler`. De belangrijkste verrijkte pathways werden weergegeven in dotplots. Ten slotte werd de TNF-signaleringsroute gevisualiseerd met het pakket `pathview` om de expressieveranderingen binnen deze pathway in kaart te brengen.

---

## Resultaten

### PCA-analyse

De PCA-analyse liet een duidelijke scheiding zien tussen de RA- en controlegroep. De eerste hoofdcomponent (PC1) verklaarde 55% van de totale variantie. Binnen beide groepen clusteren de monsters dicht bij elkaar, wat wijst op consistente genexpressieprofielen en een goede kwaliteit van de dataset zonder storende batch-effecten.

![Figuur 1: PCA Plot](figures/PCA_plot.png)

*Figuur 1. PCA-plot van de RNA-seq samples.*

### Differentiële genexpressie

De differentiële genexpressieanalyse identificeerde een groot aantal significant gereguleerde genen. In de RA-groep waren onder andere `SRGN`, `PTGFR`, `BCL2A1` en verschillende immunoglobuline-gerelateerde genen verhoogd tot expressie. Daarentegen werden `ANKRD30BL`, `MT-ND6`, `RAB3IL1` en `SLC9A3R2` significant verlaagd gevonden ten opzichte van de controlegroep.

![Figuur 2: Volcano Plot](figures/Volcano_plot.png)

*Figuur 2. Volcano plot van differentieel geëxprimeerde genen.*

### GO-verrijkingsanalyse

De GO-analyse liet zien dat voornamelijk immuungerelateerde processen verrijkt waren. Voorbeelden hiervan zijn lymfocytdifferentiatie, T-celdifferentiatie, B-celgemedieerde immuniteit en immuunrespons-gerelateerde signaalroutes. Deze resultaten suggereren dat gerichte veranderingen in de adaptieve immuunrespons een sleutelrol spelen bij RA.

![Figuur 3: GO Dotplot](figures/GO_dotplot.png)

*Figuur 3. GO Biological Process verrijkingsanalyse.*

### KEGG pathwayanalyse

De KEGG-analyse identificeerde meerdere immuun- en ontstekingsgerelateerde pathways. Belangrijke voorbeelden waren de *TNF signaling pathway*, *T cell receptor signaling pathway*, *NF-kappa B signaling pathway* en *IL-17 signaling pathway*. Deze pathways staan bekend om hun directe betrokkenheid bij ontstekingsprocessen en de pathogenese van auto-immuunziekten.

![Figuur 4: KEGG Dotplot](figures/KEGG_dotplot.png)

*Figuur 4. KEGG pathway verrijkingsanalyse.*

### TNF-signaleringsroute

De Pathview-visualisatie van de TNF-signaleringsroute liet zien dat meerdere genen binnen deze route differentieel geëxprimeerd waren. Zowel opgereguleerde als neer-gereguleerde genen waren zichtbaar binnen verschillende onderdelen van de route, met een duidelijke activatie van downstream inflammatoire cytokines. Dit ondersteunt het belang van de TNF-signaleringsroute bij RA.

![Figuur 5: TNF-signaling pathway](figures/TNF_signaling_pathway.png)

*Figuur 5. Visualisatie van differentieel geëxprimeerde genen binnen de TNF-signaleringsroute.*

---

## Conclusie

In deze studie werd met behulp van RNA-sequencing onderzocht welke genexpressieveranderingen geassocieerd zijn met reumatoïde artritis. De resultaten laten zien dat RA gepaard gaat met duidelijke verschillen in genexpressie ten opzichte van gezonde controles.

De GO- en KEGG-analyses tonen aan dat vooral immuunresponsen, lymfocytactivatie en ontstekingsgerelateerde pathways betrokken zijn bij de ziekte. Daarnaast bevestigt de verrijking van de TNF-signaleringsroute het fundamentele belang van deze pathway in de pathogenese van RA.

De sterke betrokkenheid van de TNF-signaleringsroute sluit naadloos aan bij de klinische praktijk, waarin TNF-remmers succesvol worden toegepast als biologicals bij de behandeling van RA. Deze resultaten laten zien dat transcriptomische analyse een waardevolle methode is om de moleculaire mechanismen van auto-immuunziekten beter te begrijpen.

---

##Gebruik van GitHub

Tijdens dit transcriptomicsproject heb ik onderzoeksgegevens op een gestructureerde en reproduceerbare manier beheerd met behulp van GitHub en GitHub Desktop. Om overzicht te behouden heb ik een duidelijke mappenstructuur gebruikt, bestaande uit de mappen data, scripts, results en figures. Hierdoor zijn invoerbestanden, analysescripts, resultaten en visualisaties eenvoudig terug te vinden.

Voor versiebeheer heb ik GitHub Desktop gebruikt. Tijdens het project werden wijzigingen opgeslagen met commits, zodat eerdere versies behouden bleven en aanpassingen konden worden gevolgd. Dit maakte het mogelijk om georganiseerd te werken en veranderingen overzichtelijk te documenteren.

De oorspronkelijke dataset is niet aangepast. Verwerkte bestanden en resultaten zijn apart opgeslagen om onderscheid tussen ruwe data en analyse-output te behouden. Om de reproduceerbaarheid van het onderzoek te waarborgen zijn alle gebruikte scripts, invoerbestanden en resultaten opgenomen in de repository. Daarnaast beschrijft de README de uitgevoerde analyses en de structuur van het project.

GitHub werd gebruikt als centraal platform voor het beheren en documenteren van het project. Hierdoor zijn alle onderdelen van de analyse transparant beschikbaar. Door deze werkwijze kan een andere gebruiker de scripts uitvoeren, de resultaten controleren en de analyse opnieuw reproduceren. Deze aanpak sluit aan bij de principes van transparant en betrouwbaar wetenschappelijk onderzoek.

---

## Referenties

Conesa, A., Madrigal, P., Tarazona, S., Gomez-Cabrero, D., Cervera, A., McPherson, A., Szcześniak, M. W., Gaffney, D. J., Elo, L. L., Zhang, X., & Mortazavi, A. (2016). A survey of best practices for RNA-seq data analysis. Genome Biology, 17(1), 13. https://doi.org/10.1186/s13059-016-0881-8

Smolen, J. S., Aletaha, D., & McInnes, I. B. (2016). Rheumatoid arthritis. The Lancet, 388(10055), 2023–2038. https://doi.org/10.1016/S0140-6736(16)30173-8

Zhang, F., Wei, K., Slowikowski, K., Fonseka, C. Y., Rao, D. A., Kelly, S., Goodman, S. M., Tabechian, D., Hughes, L. B., Salomon-Escoto, K., Watts, G. F. M., Jonsson, A. H., Rangel-Moreno, J., Meednu, N., Rozo, C., Apruzzese, W., Eisenhaure, T. M., Lieb, D. J., Boyle, D. L., Mandelin, A. M., II, et al. (2019). Defining inflammatory cell states in rheumatoid arthritis joint tissues by integrating single-cell transcriptomics and mass cytometry. Nature Immunology, 20(7), 928–942. https://doi.org/10.1038/s41590-019-0378-1