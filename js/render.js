import * as param from "./const.js"

const fixDiv = document.querySelector('.app__fix');
const tbody = document.querySelector('#data');

fixDiv.innerHTML = `
    <div><b>t1</b>${param.t1} [c]</div>
    <div><b>t2</b>${param.t2} [c]</div>
    <div><b>thet_torch</b>${param.thet_torch} [рад]</div>
    <div><b>thet_2</b>${param.thet_2} [рад]</div>
`
const tr_iter = (t,el) => {
    return `
        <td class="app__ceil">${t.toFixed(5)}</td>
        <td class="app__ceil">${el[0].toFixed(4)}</td>
        <td class="app__ceil">${el[1].toFixed(4)}</td>
        <td class="app__ceil">${el[2].toFixed(4)}</td>
        <td class="app__ceil">${el[3].toFixed(4)}</td>
        <td class="app__ceil">${el[4].toFixed(4)}</td>
        <td class="app__ceil">${el[5].toFixed(4)}</td>
        <td class="app__ceil">${el[6].toFixed(4)}</td>
        <td class="app__ceil">${el[7].toFixed(4)}</td>
        <td class="app__ceil">${el[8].toFixed(4)}</td>
        <td class="app__ceil">${el[9].toFixed(4)}</td>
        <td class="app__ceil">${el[10].toFixed(4)}</td>
    `
}

const html = (t, el) => {
    const tr = document.createElement('tr')
    tr.insertAdjacentHTML('beforeend', tr_iter(t, el))
    tbody.insertAdjacentElement('beforeend', tr)
}

export {
    tr_iter,
    html
}