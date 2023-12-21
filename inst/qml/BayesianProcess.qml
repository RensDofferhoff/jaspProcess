//
// Copyright (C) 2023 University of Amsterdam and Netherlands eScience Center
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public
// License along with this program.  If not, see
// <http://www.gnu.org/licenses/>.
//

import QtQuick
import QtQuick.Layouts
import JASP
import JASP.Controls
import "./common" as Common

Form
{
	Common.VariablesForm {}

    Section
    {
        title: qsTr("Models")
        columns: 1

        TabView
        {
            id:				models
            name:			"processModels"
            maximumItems:	10
            newItemName:	qsTr("Model 1")
            optionKey:		"name"

            content: Group
            {
                Group {
                    anchors.left: 		parent.left
                    anchors.margins: 	jaspTheme.contentMargin
                    
					Common.InputType
					{
						id: inputType
						modelName: rowValue
					}

                    Common.Separator {}

                    Common.InputVariables {
						visible: inputType.value == "inputVariables"
						adjustedWidth: models.width - 2 * jaspTheme.contentMargin
						colWidth: (models.width - 3 * 40 * preferencesModel.uiScale) / 4
					}

                    Common.InputModelNumber {
						visible: inputType.value == "inputModelNumber"
						adjustedWidth: models.width - 2 * jaspTheme.contentMargin
					}

					Group
					{
                        id: opts
						title: 		qsTr("Options for %1").arg(rowValue)
						columns: 	3
                        
						Group
						{
                            title: qsTr("Residual Covariances")

							CheckBox
							{
								name: "independentCovariances"
								label: qsTr("Independent variables")
								checked: independentCovariancesForAllModels.checked
							}
							CheckBox
							{
								name: "mediatorCovariances"
								label: qsTr("Mediators")
								checked: mediatorCovariancesForAllModels.checked
							}
						}

						Group
						{
                            title: qsTr("Parameter Estimates")
							// columns: 	1
							CheckBox
							{
								name: "pathCoefficients"
								label: qsTr("Paths")
								checked: pathCoefficientsForAllModels.checked
							}
							CheckBox
							{
								name: "mediationEffects"
								label: qsTr("Mediation")
								checked: mediationEffectsForAllModels.checked
							}
							CheckBox
							{
								name: "totalEffects"
								label: qsTr("Total")
								checked: totalEffectsForAllModels.checked
							}
							CheckBox
							{
								name: "residualCovariances"
								label: qsTr("Residual covariances")
								checked: residualCovariancesForAllModels.checked
							}
						}

						Group
						{
                            title: qsTr("Path Plots")
							columns: 	1
							CheckBox
							{
								name: "conceptualPathPlot"
								label: qsTr("Conceptual")
								checked: conceptualPathPlotsForAllModels.checked
							}
							CheckBox
							{
								name: "statisticalPathPlot"
								label: qsTr("Statistical")
								checked: statisticalPathPlotsForAllModels.checked
							}
						}
					}
                }
            }
        }
    }

    Section
    {
        title: qsTr("Options")
        columns: 3

		Group
		{
			IntegerField
			{
				name:			"mcmcBurnin"
				id:				warmup
				label:			qsTr("Burnin")
				defaultValue:	500
				min:			100
			}

			IntegerField
			{
				name:			"mcmcSamples"
				label:			qsTr("Samples")
				defaultValue:	1000
				min:			parseInt(warmup.value) + 100
			}

			IntegerField
			{
				name:			"mcmcChains"
				label:			qsTr("Chains")
				defaultValue:	3
				min:			1
			}

			SetSeed {}
		}
		Group
		{
			CheckBox 
			{
				label: qsTr("Parameter labels")
				name: "parameterLabels"
			}
			CheckBox 
			{
				label: qsTr("Lavaan syntax")
				name: "syntax"
			}
		}

		Common.StandardizedEstimates {}

		CIField 
		{
			text: qsTr("Credible intervals")
			name: "ciLevel"
		}
	}

	Section
	{
		title:		qsTr("MCMC diagnostics")

		VariablesForm
		{
			preferredHeight: 200 * preferencesModel.uiScale

			AvailableVariablesList
			{
				name:	"mcmcDiagnosticsAvailableTerms"
				title:	qsTr("Model terms")
				source:	[ "covariates", "factors" ]
			}

			AssignedVariablesList
			{
				singleVariable:	true
				name:			"mcmcDiagnosticsHorizontal"
				title:			mcmcDiagnosticsType.currentValue == "scatterplot" ? qsTr("Horizontal axis") : qsTr("Plotted term")
			}

			AssignedVariablesList
			{
				singleVariable:	true
				name:			"mcmcDiagnosticsVertical"
				title:			qsTr("Vertical axis")
				visible:		active
				
				property bool active:	mcmcDiagnosticsType.currentValue == "scatterplot"
				onActiveChanged:		if (!active && count > 0) itemDoubleClicked(0)
			}
		}

		DropDown
		{
			name:	"mcmcDiagnosticsType"
			id:		mcmcDiagnosticsType
			label:	qsTr("Plot type")
			values:
			[
				{ label: qsTr("Traceplot"),			value: "traceplot"},
				{ label: qsTr("Scatterplot"),		value: "scatterplot"},
				{ label: qsTr("Histogram"),			value: "histogram"},
				{ label: qsTr("Density"),			value: "density"},
				{ label: qsTr("Autocorrelations"),	value: "autocorrelation"}
			]
		}
	}

	Common.PlotOptions {}

	Section 
    {
        id: advanced
        title: qsTr("Advanced")
        columns: 1

        Group
        {
            title: qsTr("Set for All Models")
            columns: 4

			Group
			{
                title: qsTr("Residual Covariances")
				CheckBox
                {
                    id:			independentCovariancesForAllModels
                    name: 		"independentCovariancesForAllModels"
                    label: 		qsTr("Independent variables")
                }
				CheckBox
                {
                    id:			mediatorCovariancesForAllModels
                    name: 		"mediatorCovariancesForAllModels"
                    label: 		qsTr("Mediators")
                }
			}

            Group
            {
                title: qsTr("Parameter Estimates")
				columns: 1

                CheckBox
                {
                    id:			pathCoefficientsForAllModels
                    name: 		"pathCoefficientsForAllModels"
                    label: 		qsTr("Paths")
                    checked: 	true
                }
                CheckBox
                {
                    id:			mediationEffectsForAllModels
                    name: 		"mediationEffectsForAllModels"
                    label: 		qsTr("Mediation")
                    checked: 	true
                }
				CheckBox
                {
                    id:			totalEffectsForAllModels
                    name: 		"totalEffectsForAllModels"
                    label: 		qsTr("Total")
                    checked: 	true
                }
				CheckBox
                {
                    id:			residualCovariancesForAllModels
                    name: 		"residualCovariancesForAllModels"
                    label: 		qsTr("Residual covariances")
                }
            }

            Group
            {
                title: qsTr("Path Plots")
                CheckBox
                {
                    id:			conceptualPathPlotsForAllModels
                    name: 		"conceptualPathPlotsForAllModels"
                    label: 		qsTr("Conceptual")
                    checked: 	true
                }
                CheckBox
                {
                    id:			statisticalPathPlotsForAllModels
                    name: 		"statisticalPathPlotsForAllModels"
                    label: 		qsTr("Statistical")
                }
            }
		}

		Common.ModerationProbes {}
    }
}
